module PacBio.Hmm.Recursions

open System
open System.Collections.Generic 
open PacBio.Utils
open System.Diagnostics
open PacBio.Align
open PacBio.Hmm.Utils
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics
open System.Linq

(* PacBio HMM Overview

  Important types:

   ednaParams -- high level 'physical' parameters we want to estimate
   transEmDists -- compact representation of transistion and emission distributions of HMM states
   int model -- an   HMM model, including a template sequence, that generates integer observation (observation
                are 1,2,3,4 for the 4 colour channels, with 0 resevered for the null observation
   sufficientStatistics -- counts of the observation of different event types.


  The HMM equations for transition-tied observations with nulls:
  HMM Forward variable:
  q_i -> State sequence, indexed by time variable i
  S_i -> List of available states, indexed by i, 
  S_N -> final state
 
  Transition probabilities:
  a_{ij} = P(q_t = S_i | q_{t-1} = S_j, L)
  Observation probabilities:
  b_{ij}(O) = P(O | q_t = S_i, q_{t-1} = S_j, L)
 
  Forward variable:
  A_{n,t}(i) = P(O_1, O_2, ... O_n, q_t = S_i | L) 
 
  Induction of forward variable:
  P(O_1 O_2 .. O_n, q_t = S_i | L) =    \Sum_j [ P(O_1, O_2, .. O_{n-1}, q_{t-1} = S_j | L) a_{ij} b_{ij}(O_n) ] 
                                      + \Sum_j [ P(O_1, O_2, .. O_{n}, q_t = S_j | L) a_{ij} b_{ij}(null) ] 
  
  A_{n,t}(i) = \Sum_j [ A_{n-1, t-1}(j) a_{ij} b_{ij}(O_n) ] + \Sum_j [ A_{n, t-1} a_{ij} b_{ij}(null) ]
 
  Probability of Seq: P(O_1, O_2, ... O_n, q_t = S_N)
 
  Viterbi D variable: Probability of most likely state sequence
 
  D_{n,t}(i) = max_{q_1 .. q_t}  P(q_1 .. q_t = i, O_1 .. O_n | L)
                 
  No null observations:
             = max_j ( D_{n-1, t-1}(j) a_{i,j} b_{i,j}(O_n) )
 
  With null observations:
             = max {  max_j( D_{n, t-1}(j) a_{i,j} b_{i,j}(null) ), max_j( D_{n-1, t-1}(j) a_{i,j} b_{i,j}(O_n) )  } 
 *)




// You either see a numbered base or nothing     
type 'a observation =
    | Obs of 'a
    | Null   
     
// States are numbered and may have a list of transitions and their probabilities.
// a 'Final' state terminates observations
type 'a state = 
    | State of (int * dist<'a transition>) 
    | Junk of (int * dist<'a transition>)
    | Final of int
and 'a transition = 
    Transition of (int * int * dist<'a observation>)

//let (|ToVect|) (a:'a[]) = Vector.Generic.ofArray a
//let (|ToArr|)  (a:Vector<'a>) = Vector.Generic.toArray a


// final 'high-level' Edna parameters -- this is what is output
// see helper functions associated with this type below
type ednaParams = { 
    numbases:int; 
    insert:float[,]; 
    miscall:float[,]; 
    dark:float[]; 
    merge:float[] }

type ednaParams with
    // Methods for blitting this type to and from linear arrays.  Comes in handy when trying to interface with 
    // optimizers etc.
    static member toArray (e:ednaParams) = 
        (Array2D.toArray e.miscall) +> (Array2D.toArray e.insert) +> e.dark +> e.merge |> Seq.toArray
    
    static member arrayColumnNames() =
        let mc = Array2D.init 4 4 (fun target src -> "Miscall" + src.ToString() + "To" + target.ToString())
        let ins = Array2D.init 4 4 (fun target src-> "InsertOn" + src.ToString() + "Of" + target.ToString())
        let dark = Array.init 4 (fun ch -> "Dark" + ch.ToString())
        let merge = Array.init 4 (fun ch -> "Merge" + ch.ToString())

        (Array2D.toArray mc) +> (Array2D.toArray ins) +> dark +> merge |> Seq.toArray

    static member ofArray n (a : float[]) : ednaParams =
        {   numbases = n;
            miscall = Array2D.ofArray n n a.[0..n*n-1]
            insert = Array2D.ofArray n n a.[n*n..n*n+n*n-1]
            dark = a.[2*n*n..2*n*n+(n-1)]
            merge = a.[2*n*n+n..2*n*n+2*n-1]
        }
    
    static member starter =
        {
            numbases = 4;
            insert = Array2D.create 4 4 0.05 ;
            dark = Array.create 4 0.05;
            miscall = Array2D.init 4 4 (fun i j -> if i = j then 0.98 else 0.02 / 3.0)
            merge = Array.create 4 0.12
        }
        
    (* Each row defines a constraint on the flattened version of a byBaseErrorParams.
       We make one constraint for each column of the miscall matrix that forces it to add up to 1.
       The constrained solver will only ever choose parameter vectors where Ax=b. *)
    static member constraint_array2_A n =        
        let cc r c = if c >= n*r && c < n*(r+1) then 1.0 else 0.0          
        Array2D.init n (n*(2*n+2)) cc

    static member constraint_matrix_A n =        
        let cc r c = if c >= n*r && c < n*(r+1) then 1.0 else 0.0          
        Matrix.Build.Dense(n, (n*(2*n+2)), cc)

    static member constraint_array_b n = Array.create n 1.

    // FIXME - not correct
    static member total_error (e : ednaParams) =
        let nonMiscall = (e.dark |> Seq.average) + (e.merge |> Seq.average) / 4.0   
        let mc = Array.init e.numbases (fun i -> 1.0 - e.miscall.[i,i]) |> Seq.average
        nonMiscall + mc


// a compact representation of the HMM transition and emission matrices.
// the transition probabilities for each base are encoded in pStay and pMerge.
// the emission probabilites for incorporation is in moveDists and the emission probabilities 
// for insertions are in staysDists.
type transEmDists = { pStay: float[]; pMerge: float[]; moveDists: float[][]; stayDists: float[][] } with

    static member toArray (dist : transEmDists) = 
        let n = 4 + 4*5 + 4*5 + 4
        let arr = Array.create n 0.0

        let c = ref 0
        let stuff (a : float[]) = 
            for i in 0 .. a.Length - 1 do
                arr.[!c + i] <- a.[i]
            c := (!c) + a.Length

        // stuff many
        let sm (a: float[][]) = Array.iter (fun v -> stuff v) a
        
        stuff dist.pStay
        sm dist.stayDists
        sm dist.moveDists
        stuff dist.pMerge

        arr

    static member ofArray (distArr : float[]) =
        let start = ref 0

        let take n =
            let r = distArr.[!start..(!start+n-1)]
            start := !start + n
            r

        let multiTake outer inner = Array.init outer (fun i -> take inner)

        let pStay = take 4;
        let stayDists = multiTake 4 5;
        let moveDists = multiTake 4 5;
        let pMerge = take 4;

        {
            pStay = pStay
            stayDists = stayDists
            moveDists = moveDists
            pMerge = pMerge
        }

    static member has_changed absTol relTol oldpars newpars =
        let a1 = transEmDists.toArray oldpars
        let a2 = transEmDists.toArray newpars
        
        let (relIdx, maxRel) = Seq.init a1.Length (fun idx -> (idx, if a2.[idx] < absTol then 0.0 else Math.Abs((a2.[idx]-a1.[idx])/a2.[idx]))) |> Seq.maxBy snd
        let (absIdx, maxAbs) = Seq.init a1.Length (fun idx -> (idx, Math.Abs(a2.[idx]-a1.[idx]))) |> Seq.maxBy snd

        if maxAbs > absTol || maxRel > relTol then
            true
        else
            false
    
    // static member starter is defined below -- it provides a 'default' transEmDists value


// represents the transition and emission matrices and the associated error estimates
type transEmDistsErr = { dists : transEmDists; errors: transEmDists }


// An HMM model with some important metadata
type 'a model = { states: unit -> 'a state list; template:int[]; numbases:int; baseConvert: 'a->int; symbolMap: 'a->int; dists: transEmDists }

// CoCheckSumNearOneoat[] to a vector, leaving a space at the beginning for the dark observation
let PadVec (arr:float[]) = Vector<float>.Build.Dense( (arr.Length + 1), (fun i -> if i=0 then 0. else arr.[i-1]) )

let PadMat (arr:float[,]) = 
    let m = DenseMatrix.init (arr.GetLength(0)+1) (arr.GetLength(1)+1) (fun i j -> if i=0 || j=0 then 0. else arr.[i-1,j-1])
    m.[0,0] <- 1.0
    m



// convert ednaParams to transEmDists using the mapping of physical error parameters into
// HMM parameters. The reverse translation is done below.
let ednaParamsToTransEmDists (e:ednaParams) =
       
    let _insertMat = PadMat e.insert        
    let _miscallMat = PadMat e.miscall
   

    // Miscall matrix has to be stochastic, insert matrix doesn't
    let miscallMat = fixup _miscallMat
    let insertMat = _insertMat
    
    let dark = PadVec e.dark
    let merge = PadVec e.merge

    let p2e p = p / (1. - p)
    let e2p e = e / (1. + e)
    
    let insertMatExp = Matrix.map p2e insertMat
    let channels = [|1..e.numbases|]
    let obs = [|0..e.numbases|]

    // Make normalized branch and stick probabilities
    // Be careful to actually normalize these things properly

    let eStay = channels |> Array.map (fun cog -> [| for o in obs -> insertMatExp.[o,cog] |] |> Array.sum)
    let pStay = eStay |> Array.map e2p    
    let vStay = channels |> Array.map (fun cog -> obs |> Array.map (fun o -> insertMatExp.[o,cog] / eStay.[cog-1]))

    // Merge probability defined as a fraction of 2-mer incorporations - so scale down by 1-pStay
    let pMerge = channels |> Array.map (fun i -> (1.0 - pStay.[i-1]) * merge.[i])

    // Fill out the dark matrix that takes an observation matrix and makes some of it go dark, according to the dark vector
    let darkMat = DenseMatrix.zero<float> (e.numbases+1) (e.numbases+1)
    Vector.toArray dark |> Array.iteri (fun idx darkVal -> darkMat.[0,idx] <- darkVal; darkMat.[idx,idx] <- 1. - darkVal)
    
    // Vector with a one in given index and zeros elsewhere
    let cog_vec bs = DenseVector.init (e.numbases+1) (fun i -> if  i=bs then 1. else 0.)
    // Constant vector
    let vconst c = DenseVector.create (e.numbases+1) c

    let moveDist bs = (miscallMat * darkMat * (cog_vec bs)).ToArray() 
    let moveDists = channels |> Array.map moveDist
    
    let stayDist bs = vStay.[bs-1].ToArray()
    let stayDists = channels |> Array.map stayDist

    // when merging, we must emit the cognate base
    let mergeDists = channels |> Array.map (fun bs -> (cog_vec bs).Clone())

    // Check that the error probabilities are valid
    let cp p = if p >= 0. && p < 1. then true else false
    // Check that all elements of a vector have valid probabilites
    let checkVec = Vector.forall cp  
    // Check that a matrix preserves normalized probabilities         
    if not (isStochastic 0.001 miscallMat) then failwith("Miscall matrix must be stochastic (column sums == 1)")
    if not (isStochastic 0.001 darkMat) then failwith("Problem constructing dark matrix")
         
    { pStay = pStay; pMerge = pMerge; stayDists = stayDists; moveDists =  moveDists }


type transEmDists with
    static member starter = ednaParamsToTransEmDists ednaParams.starter

    
let MinP = Math.Sqrt(Double.Epsilon)
let MaxP = 1.0 - 0.0000001


// Convert the transEmDists representation to the ednaParams representation
let transEmDistsToEdnaParams (dists : transEmDistsErr) =
    let numbases = dists.dists.pStay.Length

    let distArray = transEmDists.toArray dists.dists
    let errorArray = transEmDists.toArray dists.errors

    if distArray |> Array.exists Double.IsNaN then 
        printfn "Got NaN in Param convert"

    let startGenParamsVect = DenseVector.raw (ednaParams.toArray ednaParams.starter )

    let fmin (trialGenParamsVect : Vector<float>) = 
        let  trialGenParamsArray =  trialGenParamsVect.Map((fun v -> Math.Max(MinP, Math.Min(MaxP, v)))).ToArray()
        // lm may not respect the lb and ub that you give it, so fix the argument.
        let trialGenParams = ednaParams.ofArray numbases trialGenParamsArray
        let trialDistArray = transEmDists.toArray (ednaParamsToTransEmDists trialGenParams)
        let diffVect = DenseVector.init distArray.Length (fun i -> (distArray.[i] - trialDistArray.[i]) / errorArray.[i])
        diffVect

    // Use a constrained LM solver to enforce a normalized miscall matrix
    let A = ednaParams.constraint_matrix_A numbases
    let b = ednaParams.constraint_array_b numbases
    let lb = DenseVector.create startGenParamsVect.Count MinP
    let ub = DenseVector.create startGenParamsVect.Count MaxP
        
    // lmsove_con does not respect the lb and ub that you give it, so fix the argument.
    let eps = 1e-7
    let opts = { lm2.OptimizationOptions.standard(50) with Eps1 = eps; Eps2 = eps; Eps3 = eps }
    let ra = lm2.lmBoxConstrains (Some(lb)) (Some(ub)) A (Vector.ofArray b) fmin startGenParamsVect opts`

    //printfn "solved: %A" ra

    let flatednaParams = ra.p |> Vector.map (fun v -> Math.Max(MinP, Math.Min(MaxP, v))) |> Vector.toArray
    ednaParams.ofArray numbases flatednaParams

   
// Build a simple left-to-right HMM for sequencing with a given 
// stick rate, branching fraction, 'missed observation', and miscall fraction
// The input vectors (stick, branch, and dark) are interpreted as per-base probabilities.  
// The miscall matrix is interpreted as a stochastic mixing matrix over the n nucleotides.
let makeModelFromTransEmDists (nbases : int) (d : transEmDists) (template : int array) =
    
    let _makeStates() =
        let nf = float(nbases)
        let sequence = List.ofArray template
        let m = List.length sequence
    
        // Make sure all the elements of the sequence are in the specified range
        //if Array.forall (ismember {1 .. nbases}) template then () else failwith("Sequence elements must be in range 1 to n")
        for i = 0 to (template.Length - 1) do
            if template.[i] < 1 || template.[i] > nbases then failwith("Sequence elements must be in range 1 to n")
        done
    
        // Convert distribution into a list of (probability:float, base:obserservation) tuples
        let mkObservationList observation_dist = 
            observation_dist |> List.ofArray |> 
            List.mapi (fun b v -> if b=0 then (v, Null) else (v, Obs(b)))

        // We apply the dark matrix first -- this is somewhat arbitrary, but dark dyes should remain dark, and
        // not get a new lease on life when the are miscalled to a brighter channel for example
        let _stayObs bs = d.stayDists.[bs] |> mkObservationList  
        let stayObs = memoize _stayObs
    
        // The observation vector if we move
        let _moveObs bs = d.moveDists.[bs] |> mkObservationList
        let moveObs = memoize _moveObs

        //let _mergeObs bs = d.mergeDists.[bs] |> mkObservationList
        let _mergeObs bs = [0..4] |> List.mapi (fun b v -> if b = 0 then (0.0, Null) else (if bs + 1 = b then (1.0, Obs(b)) else (0.0, Obs(b))))
        let mergeObs = memoize _mergeObs


        //  A state consists of a number, and a list of (probability, state) tuples
        let emitState n bs mergeable = 
            let sObs = stayObs (bs-1)
            let mObs = moveObs (bs-1)
        
            let ps = d.pStay.[bs-1]               

            if mergeable then
                let mgObs = mergeObs (bs - 1)
                let pm = (1.0-ps) * d.pMerge.[bs-1]

                State(n, [ (ps, Transition(n, n, sObs)) ; ( (1. - ps - pm), Transition(n, n+1, mObs)); (pm, Transition(n, n+2, mgObs)) ])
            else
                State(n, [ (ps, Transition(n, n, sObs)) ; ( (1. - ps), Transition(n, n+1, mObs)) ])
    
        // Generate a list of states based on the sequence
        let seqLength = template.Length

        let doEmit n bs =
            if n < seqLength - 1 && sequence.[n+1] = bs then
                emitState n bs true
            else
                emitState n bs false

        let model = (template |> Seq.mapi (fun n bs ->  doEmit n bs) |> Seq.toList) @ [Final(List.length sequence)]
    
        // Do some sanity checks on the model please!!
        if verify_model model then () else Console.WriteLine (sprintf "%A" model); failwith("Couldn't make a proper model")
        model

    {states = _makeStates; template = template; numbases = nbases; baseConvert = id; symbolMap = id; dists = d} : int model


// Hepler function to build an HMM model.  Conver the edna params to transEmDists, then make model.
let makeModelFromEdnaParams tpl (ednaParams : ednaParams) = 
    makeModelFromTransEmDists ednaParams.numbases (ednaParamsToTransEmDists ednaParams) tpl
    

// Generate sample observation sequences from the given HMM
let sampleModel (m : 'a model) =
    let model = m.states() |> Seq.toArray

    let getObsandNext transes = 
        match (sampleDist transes) with Transition(source, dest, obs) -> (sampleDist obs, (model.[dest]))

    let rec makeSequence state seq = 
        match state with
            | Final(_) -> seq
            | Junk(_) -> seq
            | State(_, transes) -> 
                match getObsandNext transes with
                    | (Null, next) -> (makeSequence next seq)
                    | (Obs(b), next) -> (makeSequence next (b :: seq))
            
    List.toArray (List.rev (makeSequence (model.[0]) []))
  


type hmmEvaluators(_p : LFloat, _alpha : IBandedMatrix<LFloat>, _beta : IBandedMatrix<LFloat>, _viterbi : unit -> AlignCell[], dispose: unit -> unit) =
    interface IDisposable with
        member x.Dispose() = dispose()
    member x.p = _p
    member x.alpha = _alpha
    member x.beta = _beta
    member x.viterbi = _viterbi


let toFloatVector (arr : float[]) =
    let fa = new ConsensusCore.FloatVector(arr.Length)
    for i = 0 to arr.Length - 1 do
        fa.Add(float32 arr.[i])
    done

    fa
    
let toFloatVector2 (arr : float[][]) =
    let fa = new ConsensusCore.FloatVector(arr |> Seq.sumBy (fun a -> a.Length))

    let p = ref 0

    for i = 0 to arr.Length - 1 do
        let a = arr.[i]
        for j = 0 to a.Length - 1 do

            fa.Add(float32 a.[j])
            incr p
        done
    done

    fa

let toIntVector (arr : int[]) =
    let fa = new ConsensusCore.IntVector(arr.Length)
    for i = 0 to arr.Length - 1 do
        fa.Add(arr.[i])
    done

    fa

let distsToEdnaParams (dists : transEmDists) =
    let arrays = ref [] 

    let cnv (a: float[]) =
        let fa = (toFloatVector a)
        arrays := fa :: !arrays
        fa

    let cnv2 (a : float[][]) =
        let fa = (toFloatVector2 a)
        arrays := fa :: !arrays
        fa
         
    let p = 
        new ConsensusCore.EdnaModelParams(
            cnv dists.pStay,
            cnv dists.pMerge,
            cnv2 dists.moveDists,
            cnv2 dists.stayDists)

    let disp = 
        { new IDisposable with 
            member this.Dispose() = !arrays |> List.iter (fun i -> i.Dispose()) }
    (p, disp)


let scoreDiff = ref 35.0f

 // Compute the probability of the sequence given the model
let setupHmmEvaluators (model : 'a model) (read : int[]) =
    let M = model.template.Length + 1
    let N = read.Length   

    let dists = model.dists
    let nb = model.numbases

    let tpl = model.template |> Array.map byte
    
    // evaluator parameters
    let (pars, parsDispose) = distsToEdnaParams model.dists
    
    // Read data
    let stringRead = new String(read |> Array.map char)
    use intRead = (toIntVector read)
    use features = new ConsensusCore.ChannelSequenceFeatures(stringRead, intRead)

    // template and evaluator
    let stringTpl =  new String(tpl |> Array.map char)    
    use intTpl = (toIntVector (tpl |> Array.map int))
    use eval = new ConsensusCore.EdnaEvaluator(features, stringTpl, intTpl, pars)
    
    //validateEval features stringTpl intTpl pars

    // banding options    
    use bandOptions = new ConsensusCore.BandingOptions(4, !scoreDiff)

    // put the recursor and the scorer together
    use recursor = new ConsensusCore.SparseSseEdnaRecursor(15,bandOptions)
    let scorer = new ConsensusCore.SparseSseEdnaMutationScorer(eval, recursor)

    let alpha = scorer.Alpha()
    let beta = scorer.Beta()
    let p0 = alpha.Get(N,(M-1)) |> LFloat.Explicit

    let makeQuick (mat : ConsensusCore.SparseMatrix) = 
        { new IBandedMatrix<LFloat> with
             member x.Item
                with get(i,j) = LFloat.Explicit(mat.Get(i,j))
                and set(i,j) (v:LFloat) = () // mat.[i,j] <- v.Value
             member x.Rows with get() = mat.Rows()
             member x.Cols with get() = mat.Columns()
             member x.RowRangeForCol(i) = new IndexRange(0,N) //mat.UsedRowRange(i)
             member x.ColRangeForRow(i) = new IndexRange(0,M) //mat.ColRangeForRow(i)
             member x.ToArray() = failwith("not implemented") }

    let viterbiCell () = [||]
        //let vc = Viterbi.ViterbiRecursion(tpl.ToArray(), passModel)
        //vc.ToArray()

    let dispose() = 
        scorer.Dispose()

    new hmmEvaluators(p0, makeQuick alpha, makeQuick beta, viterbiCell, dispose)


type tempateMutation = Substitution of (int * int) | Deletion of int | Insertion of int * int

// Utilities for adding LFloat
let ldsum (a : #seq<LFloat>) = Seq.fold (fun (a : LFloat) b -> a + b) LFloat.Zero a
let ldsuma = Array.fold (fun (a : LFloat) b -> a + b) LFloat.Zero
