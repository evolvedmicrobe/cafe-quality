module PacBio.Hmm.ExpectationMaximization

open System
open System.Collections.Generic 
open Microsoft.FSharp.Math
open PacBio.Utils
open PacBio.Hmm.Recursions    
open Microsoft.FSharp.Math         
open PacBio.Align

open System.Diagnostics

let ptwise_combine f (a: 'a[]) (b: 'b[]) =
    if a.Length <> b.Length then failwith("Arrays dimension mismatch")
    Array.init (a.Length) (fun i -> f a.[i] b.[i])

let arraylsum n (arrays: #seq<LFloat[]>) = 
    let array_cache = arrays |> Seq.cache
    let res = Array.create n (LFloat.Zero)
    for a in array_cache do
        for i in 0 .. res.Length - 1 do
            res.[i] <- res.[i] + a.[i]
        done
    done
    res
    
let arraylsuma n (arrays: LFloat[] array) = 
    let res = Array.create n (LFloat.Zero)
    for a in arrays do
        for i in 0 .. res.Length - 1 do
            res.[i] <- res.[i] + a.[i]
        done
    done
    res    
    
    
let arraysub a b = ptwise_combine (-) a b
let arrayadd a b = ptwise_combine (-) a b


// The sufficient statistics counts the expected number of events of different types when analyzing a read and template with a given HMM model 
type sufficientStatistics = 
    { 
        eta_stay:LFloat[];      // Number of stay transitions per base
        eta_move:LFloat[];      // Number of move transitions per base
        eta_no_merge:LFloat[];  // Number of move transitions per base in mergable positions
        eta_merge:LFloat[];     // Number of merge transitions per base in mergable positions
        etao_stay:LFloat[][];   // Counts of observations when in a stay loop, per base
        etao_move:LFloat[][];   // Counts of observations when in a move trans, per base
        etao_merge:LFloat[][]   // Counts of observations when in a merge trans -- merges are currently constrained to always emit the cognate signal
    }

type sufficientStatistics with
    // Sufficient statistics add simply.
    static member ( + ) (a, b) = 
        let n = a.eta_stay.Length
        {
            eta_stay     = arraylsum n [a.eta_stay; b.eta_stay];
            eta_move     = arraylsum n [a.eta_move; b.eta_move];
            eta_no_merge = arraylsum n [a.eta_no_merge; b.eta_no_merge];
            eta_merge    = arraylsum n [a.eta_merge; b.eta_merge];
            etao_stay    = (Seq.zip a.etao_stay b.etao_stay)   |> Seq.map (fun (a,b) -> arraylsum (n+1) [a;b]) |> Seq.toArray
            etao_move    = (Seq.zip a.etao_move b.etao_move)   |> Seq.map (fun (a,b) -> arraylsum (n+1) [a;b]) |> Seq.toArray
            etao_merge   = (Seq.zip a.etao_merge b.etao_merge) |> Seq.map (fun (a,b) -> arraylsum (n+1) [a;b]) |> Seq.toArray
        }
    static member zero (numbases : int) =
        {
            eta_stay     = Array.init numbases (fun _ -> LFloat.Zero);
            eta_move     = Array.init numbases (fun _ -> LFloat.Zero);
            eta_no_merge = Array.init numbases (fun _ -> LFloat.Zero);
            eta_merge    = Array.init numbases (fun _ -> LFloat.Zero);
            etao_stay    = Array.init numbases (fun _ -> Array.init (numbases + 1) (fun _ -> LFloat.Zero));
            etao_move    = Array.init numbases (fun _ -> Array.init (numbases + 1) (fun _ -> LFloat.Zero));
            etao_merge   = Array.init numbases (fun _ -> Array.init (numbases + 1) (fun _ -> LFloat.Zero));
        }


let eps = 0
let sum s = s |> Seq.reduce (+)

type hmmPulseTag = { insert: float; miscall:float; incorporation:float; predark:float; postdark:float; merge:float }


// Test of the Edna Evaluator -- make sure it's giving normalized probabilities
let validateEval features stringTpl (intTpl : ConsensusCore.IntVector) pars =
    // Go to arbitrary tpl, read positions
    let read = [| 1;2;3;4 |]

    // Read data
    let stringRead = new String(read |> Array.map char)
    use intRead = (toIntVector read)
    use features = new ConsensusCore.ChannelSequenceFeatures(stringRead, intRead)

    // template and evaluator
    use eval = new ConsensusCore.EdnaEvaluator(features, stringTpl, intTpl, pars)

    // For all template locations all the outbound moves over all the read bases should sum to one
    let broke = ref false

    for t in 0 .. stringTpl.Length - 1 do
        let acc = ref 0.0
        let add v = acc := !acc + Math.Exp(float v)
        
        add (eval.Del(1, t))
        
        for r in 0 .. 3 do     
            add (eval.Inc(r,t))
            add (eval.Extra(r,t))
            add (eval.Merge(r,t))
        done
        
        let err = (!acc - 1.0)
        if (err > 0.1) && Debugger.IsAttached && not !broke then
            broke := true
            Debugger.Break()
    done
    

// Extract the sufficientStatistics for a read generated by a model
let extractSufficientStatisticsFast model (read : int[]) = 
    let M = model.template.Length + 1
    let N = read.Length   

    let dists = model.dists
    let nb = model.numbases

    let tpl = model.template |> Array.map byte
    
    // evaluator parameters
    let (pars, parsDispose) = distsToEdnaParams model.dists
    
    // Read data
    let stringRead = new String(read |> Array.map char) // this is junk
    use intRead = (toIntVector read)
    use features = new ConsensusCore.ChannelSequenceFeatures(stringRead, intRead)

    // template and evaluator
    let stringTpl =  new String(tpl |> Array.map char)  // this is junk
    use intTpl = (toIntVector (tpl |> Array.map int))
    use eval = new ConsensusCore.EdnaEvaluator(features, stringTpl, intTpl, pars)
    
    // Optionally double check that the EdnaEvaluator is behaving well
    //validateEval features stringTpl intTpl pars

    // banding options
    use bandOptions = new ConsensusCore.BandingOptions(4, !PacBio.Hmm.Recursions.scoreDiff)

    // put the recursor and the scorer together
    use recursor = new ConsensusCore.SparseSseEdnaRecursor(15,bandOptions)
    use scorer = new ConsensusCore.SparseSseEdnaMutationScorer(eval, recursor)

    use hmmCounter = new ConsensusCore.EdnaCounts()
  
    let po = LFloat.Explicit(scorer.Alpha().Get(N, (M-1)))

    use etaObsResult = new ConsensusCore.FloatArray(5)
    let pEtaObsResult = etaObsResult.cast()
    let channelFeature = features.Channel

    let etaObsAll i j =
        hmmCounter.DoCount(channelFeature, eval, scorer, i, j, pEtaObsResult)
        Array.init 5 (fun i -> LFloat.Explicit(etaObsResult.getitem(i)))

    let hasMerge state = 
        match state with 
            | State(n, dist) -> dist |> List.exists (fun (_, Transition(frs, tos, _)) -> frs + 2 = tos)
            | Junk(_) -> false
            | Final(_) -> false

    let hasMerge2 (template : int[]) idx =
        idx < template.Length - 1 && template.[idx] = template.[idx + 1]

    let sl = Array.init (model.numbases) (fun i -> [])
    let msl = Array.init (model.numbases) (fun i -> [])

    let idx = ref 0
    for s in 0 .. model.template.Length - 1 do
        let b = model.template.[!idx] - 1
        sl.[b] <- !idx :: sl.[b]
        
        //if hasMerge s then
        if hasMerge2 model.template (!idx) then
            msl.[b] <- !idx :: msl.[b]

        incr idx

    let stateLists = sl |> Array.map List.toArray
    let mergeStateLists = msl |> Array.map List.toArray

    let obsRange = {0..model.numbases} |> Seq.toArray    
    
    // Counts of all observation types in stay transitions, by base
    let oStay = [| for l in stateLists -> l |> Array.map (fun s -> etaObsAll s s) |> arraylsuma (nb+1) |> Array.map (fun v -> v / po) |]
    
    // Counts of all observation types in move transitions, by base
    let oMove = [| for l in stateLists -> l |> Array.map (fun s -> etaObsAll s (s+1)) |> arraylsuma (nb+1) |> Array.map (fun v -> v / po) |]
    
    // Counts of all observation types in move transitions, by base, in mergeable positions
    let oNoMerge = [| for l in mergeStateLists -> l |> Array.map (fun s -> etaObsAll s (s+1)) |> arraylsuma (nb + 1) |> Array.map (fun v -> v / po) |]

    // Counts of all observation types in merge transitions, by base, in mergeable positions
    let oMerge = [| for l in mergeStateLists -> l |> Array.map (fun s -> etaObsAll s (s+2)) |> arraylsuma (nb + 1) |> Array.map (fun v -> v / po) |]

    let aStay = oStay |> Array.map ldsuma
    let aMove = oMove |> Array.map ldsuma
    let aMerge = oMerge |> Array.map ldsuma
    let aNoMerge = oNoMerge |> Array.map ldsuma
            
    let ss =
        {eta_stay = aStay;
         eta_move = aMove;
         eta_merge = aMerge;
         eta_no_merge = aNoMerge;
         etao_stay = oStay;
         etao_move = oMove;
         etao_merge = oMerge
         }
    
    parsDispose.Dispose()

    ss
    
let findNaN (v : LFloat[][]) = Array.exists (Array.exists (float >> Double.IsNaN)) v

let extractSufficientStatistics (model : int model) (read : int[]) = 

    let M = float model.template.Length

    // Extract suff. stats
    try
        let ss = extractSufficientStatisticsFast model read

        // Verify that we are counting the correct number of forward moving bases
        let mNormal = (ss.eta_move |> ldsum |> float)
        let mMerge =  (ss.eta_merge |> ldsum |> float)

        let m = mNormal + 2.0*mMerge
        let err = Math.Abs(m-M) / (float M)

        if (Double.IsInfinity(m) || Double.IsNaN(m) || findNaN ss.etao_move || findNaN ss.etao_stay || findNaN ss.etao_merge) then
            printfn "Got NaN in suff stats: %A" ss

        // It appears that the main reason we will hit this is if there is some error in the final alpha/beta values
        // A bit of error here doesn't matter too much because when we form the transEmDists, we divide suffStats be each other
        // which should cancel out any error in 'po'
        if (err > 0.05) then
            printfn "Incorrect count of forward moves: Expected %f, Obs: %f" M m

        Some(ss)
    with
        | ex ->
            printfn "Got '%s', skipping read" ex.Message
            None

let extractSufficientStatisticsSum model reads =
    reads
    |> Seq.map (extractSufficientStatistics model)
    |> Seq.filter Option.isSome
    |> Seq.map Option.get
    |> Seq.fold (+) (sufficientStatistics.zero model.numbases)

let distsFromStats (stats:sufficientStatistics) =
    let numbases = stats.eta_stay.Length
    let templateBases = {0..numbases-1}
    
    let dm (v : LFloat) = if Single.IsNaN(v.Value) then 1e-4 else Math.Max(1e-4, double(v))

    let normalize (a : LFloat[]) =
        let a0 = a |> Array.map double
        let z0 = a0 |> Array.sum
        if z0 = 0.0 then 
            Array.create (a.Length) (1.0 / (float a.Length))
        else
            let an = a0 |> Array.map (fun v -> Math.Max(v/z0, 0.0001))
            let z = an |> Array.sum        
            an |> Array.map (fun v -> v / z)    
    
    // Estimate the stay transition probability per base
    let aii = 
        templateBases |>
        Seq.map (fun b -> dm(stats.eta_stay.[b] / (stats.eta_stay.[b] + stats.eta_move.[b]) )) |>
        Seq.toArray

    // Estimate the merge probability per hp base
    let aii2 = 
        templateBases |>
        Seq.map (fun b -> dm(stats.eta_merge.[b] / (stats.eta_merge.[b] + stats.eta_no_merge.[b]) )) |>
        Seq.toArray

    // For each starting base get the distribution of observations when going forward
    let mdist = templateBases |> Seq.map (fun b -> stats.etao_move.[b] |> normalize) |> Seq.toArray

    // Set dark merges to zero
    for t in templateBases do
        stats.etao_merge.[t].[0] <- LFloat.Zero
    done

    // Set dark inserts to zero, since we have no way of estimating them
    for t in templateBases do
        stats.etao_stay.[t].[0] <- LFloat.Zero        
    done

    let sdist = templateBases |> Seq.map (fun b -> stats.etao_stay.[b] |> normalize) |> Seq.toArray
        
    { pStay = aii; pMerge = aii2; stayDists = sdist; moveDists = mdist }


let distsErrorsFromStats (stats:sufficientStatistics) =
    let numbases = stats.eta_stay.Length
    let templateBases = {0..numbases-1}
    let obsBases = {1..numbases}
    let bases = Seq.zip templateBases obsBases
    let skipObsBase b = obsBases |> Seq.filter (fun i -> i <> b)
    
    //  Binomial Error -- take a conservative estimate -- imagine you had at least one observation of true and false
    let binErr p n0 = 
        let n = Math.Max(n0, 0.1)
        Math.Sqrt(Math.Max(p,1.0/n)*Math.Max(1.0-p, 1.0/n) / n)
    
    let d = distsFromStats stats
    let aii = d.pStay
    let mdist = d.moveDists
    let sdist = d.stayDists

    let w_aii =
        templateBases |>
        Seq.map (fun b -> binErr aii.[b] (double(stats.eta_stay.[b] + stats.eta_move.[b]))) |>
        Seq.toArray
 
    let w_aii2 =
        templateBases |>
        Seq.map (fun b -> binErr aii.[b] (double(stats.eta_no_merge.[b] + stats.eta_merge.[b]))) |>
        Seq.toArray

    let w_mdist = 
        templateBases |>
        Seq.map (fun b -> mdist.[b] |> Array.map (fun v -> binErr v (double(stats.eta_move.[b])))) |>
        Seq.toArray

    let w_sdist = 
        templateBases |>
        Seq.map (fun b -> sdist.[b] |> Array.mapi (fun i v -> if i = 0 then 100.0*(binErr 0.5 (double(stats.eta_stay.[b]))) else binErr v (double(stats.eta_stay.[b])))) |>
        Seq.toArray

    { pStay = w_aii; pMerge = w_aii2; stayDists = w_sdist; moveDists = w_mdist }

let distsAndErrorsFromStats (stats:sufficientStatistics) = 
    { dists = distsFromStats stats; errors = distsErrorsFromStats stats }



let fourCGenEstimateParams suffStats = transEmDistsToEdnaParams (distsAndErrorsFromStats suffStats)



let isFeasible arr = not (arr |> Array.exists (fun v -> v > 1.0 || v < 0.0))

let findFeasible (factor : float) (oldArr : float[]) (newArr : float[]) =

    let newTry j = Array.init newArr.Length (fun i -> newArr.[i] + (j - 1.0) * (newArr.[i] - oldArr.[i]))

    let rec jj min max =
        let mid = (max+min) / 2.0
        let midVect = newTry mid
        
        if isFeasible midVect then
            if max - mid < 0.1 then midVect else jj mid max 
        else 
            jj min mid

    let try0 = newTry factor
    if isFeasible try0 then 
        try0
    else
        jj 1.0 factor


// Use the EM algorithm to estimate the 
let estimateTransEmDists startParams absTol relTol iterations tplReadPairs  =
    let updateDistsAndErrors dists =
        let nbases = 4
        let stats (t, r) = extractSufficientStatistics (makeModelFromTransEmDists nbases dists t) r
        let all_stats = tplReadPairs |> Seq.map stats |> Seq.filter Option.isSome |> Seq.map Option.get |> Seq.fold (+) (sufficientStatistics.zero nbases)
        distsAndErrorsFromStats all_stats
    
    let rec doEstimate n (distErrs : transEmDistsErr) =
        if n = 0 then
            //Console.WriteLine(sprintf "4C EM reached max iterations: %i rounds" iterations);
            distErrs
        else
            let newDistErrs = updateDistsAndErrors distErrs.dists

            // Step doubling -- 
            let pars = (transEmDists.toArray newDistErrs.dists)
            let lastPars = (transEmDists.toArray distErrs.dists)

            let newTry = findFeasible 1.4 lastPars pars
            let jumpPars = transEmDists.ofArray newTry

            if transEmDists.has_changed absTol relTol distErrs.dists newDistErrs.dists then
                doEstimate (n-1) { newDistErrs with dists = jumpPars }
                //doEstimate (n-1) newDistErrs
            else
                //Console.WriteLine(sprintf "4C EM finished due to tolerances: %i rounds" (iterations - n + 1));
                newDistErrs
                
    let res = doEstimate iterations (updateDistsAndErrors startParams)
        
    res

let estimateTransEmDistsQuick iterations tplReadPairs = estimateTransEmDists transEmDists.starter 0.01 0.01 iterations tplReadPairs


let indmax (rv : #seq<LFloat>) = 
    let a = Seq.toArray rv
    let idx = ref 0
    let m = ref LFloat.Zero
    a |> Array.iteri (fun i v -> if v.Value > (!m).Value then idx := i; m := v) 
    !idx

let pmax = Seq.reduce max

let sw (s: int[]) = Seq.fold (fun b (a : int) -> b + a.ToString()) "" s

type LFloatMath() =
  interface INumeric<LFloat> with
    member x.Abs y = y
    member x.Add(a, b) = a + b
    member x.Multiply(a, b) = a * b
    member x.Negate a = raise (new NotImplementedException("Can't negatve an LFloat"))
    member x.Parse (str,_,_) = new LFloat(Double.Parse(str))
    member x.Sign _ = 1
    member x.Subtract(a, b) = a - b
    member x.ToString (v, str, fmt) = String.Format(fmt, str, [|(float(v) :> obj)|])
    member x.One = LFloat.One
    member x.Zero = LFloat.Zero
    member x.Compare(a,b) = if a.Value = b.Value then 0 else if a.Value > b.Value then 1 else -1
    member x.Equals(a,b) = a.Value = b.Value

GlobalAssociations.RegisterNumericAssociation (new LFloatMath())


type mutation = Insert of int * int | Deletion of int | Miscall of int * int
type lMat = Matrix<LFloat>
type lVec = Vector<LFloat>

// Returns a list of 
let templateMaximization(aii : lVec, bii : lMat, biip1 : lMat, _) (template : int[]) =
    
    let obviousChanges() =
        let acted = ref false
        let surgeryStep (aii : LFloat) (bii:#seq<LFloat>) (biip1 : RowVector<LFloat>) idx b =
            if aii.Value > (new LFloat(0.7)).Value then
                if indmax bii > 0 then
                    acted := true; [indmax bii ; b]
                else 
                    [b]
            elif pmax (biip1 |> Seq.map (fun v-> v.Value)) > biip1.[b].Value + (new LFloat(0.1)).Value then
                acted := true; 
                match indmax biip1 with 
                    | a when a = eps -> [] 
                    | newbase -> [newbase]
            else
                [b]
       
        let new_template = 
            {0 .. template.Length - 1} |> 
                Seq.map (fun i -> surgeryStep (aii.[i]) (bii.Row(i)) (biip1.Row(i)) i template.[i]) |>
                Seq.reduce (@) |> Array.ofList
        
        if (!acted) then Some(new_template) else None
    
    let killSkips() =    
        let skipcol = biip1.Column(eps)
        let skipmax = skipcol |> indmax
        let skip_reg r i = { Math.Max(0, i - r) .. Math.Min(skipcol.Length-1, i + r) }
        
        let skipRegion r = skipcol.[Math.Max(0, skipmax - r) .. Math.Min(skipcol.Length-1, skipmax + r)] |> ldsum
        if skipcol.[skipmax].Value > (new LFloat(0.5)).Value || (skipRegion 2).Value > (new LFloat(0.7)).Value then 
            Some([([ for i in {0 .. skipmax-1} -> template.[i]] @ [ for i in {skipmax+1 .. template.Length-1} -> template.[i] ]) |> Array.ofList])
        else
            None
 
    let addBranches() =
        let aiimax = aii|> indmax
        let insbase = bii.Row(aiimax) |> indmax 
        let branches = aii |> ldsum
        if branches.Value > (new LFloat(0.5)).Value && (aii.[aiimax]).Value > (new LFloat(0.1)).Value then 
            Some([([ for i in {0 .. aiimax-1} -> (template.[i]) ] @ 
                   [insbase] @ 
                   [for i in {aiimax .. template.Length - 1} -> template.[i]]) |> Array.ofList] )
        else
            None

    match obviousChanges() with
        | Some(t) -> [t]
        | None -> [killSkips; addBranches] |> Seq.choose (fun a -> a()) |> List.concat


// Store the results of evaluating a template -- a list of alternate templates, and the likelihood of this template
type emResult = { mutations : int[] list; p : LFloat }


let mMul (c : LFloat) (m : Matrix<LFloat>) = 
    let (I,J) = m.Dimensions
    for i = 0 to I-1 do
        for j = 0 to J-1 do
            m.[i,j] <- c * m.[i,j]
        done
    done
    
let vMul (c : LFloat) (m : Vector<LFloat>) = 
    let I = m.Length
    for i = 0 to I-1 do
        m.[i] <- c * m.[i]
    done