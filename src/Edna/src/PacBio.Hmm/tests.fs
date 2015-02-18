#light
module PacBio.Hmm.Test

open System
open System.Diagnostics
open System.Reflection
open System.IO
open System.Text    

open NUnit.Framework
open NUnit.Framework.Constraints
open NUnit.Framework

open Microsoft.FSharp.Math

open PacBio.Hmm.Minimize
open PacBio.Hmm.Recursions
open PacBio.Hmm.ExpectationMaximization
open PacBio.Hmm.Utils

open PacBio.Utils
open PacBio.Align

// Use a constant random seed
let rand = new Random(0)

// Some reasonable starting error parameters
let samplePars (e : float) : ednaParams =
    let em q = Array.create 4 q
    let sp v (a:float[]) = Array.init (a.Length) (fun i -> v * a.[i])

    { numbases = 4; 
        dark = sp e [|0.05; 0.03; 0.04; 0.03|]; 
        insert = Array2D.init 4 4  (fun i j -> if i=j then 0.03 * e else 0.01 * e); 
        merge = sp e [| 0.1; 0.15; 0.07; 0.1 |];  
            
        miscall = Array2D.init 4 4 (fun i j ->
            if i=j && i < 2 then 1. - 0.007 * e
            else if i=j then 1. - 0.026 * e
            else if i<2 && j<2 then 0.007 * e
            else if i>=2 && j>=2 then 0.026 * e
            else 0.)
    }


// Quick tests for utility code
[<TestFixture>]
type TestUtils() =
    
    let a = [|0;1;2;3;4;5;6;7;8;9|]

    [<Test>]
    member x.ArraySlice() =
        Assert.AreEqual(Array.Slice {0..5} a, [|0;1;2;3;4;5|])
        Assert.AreEqual(Array.Slice [1;2;3;8] a, [|1;2;3;8|])


    [<Test>]
    member x.Array2Serialize() = Assert.AreEqual(Array2D.init 2 2 (fun i j -> i+j), Array2D.ofArray 2 2 [|0;1;1;2|] ) 
    
    [<Test>]
    [<ExpectedException("System.IndexOutOfRangeException")>]
    member x.ArrayBadSlice() =
            Array.Slice [-1;3;5] a |> ignore
            
    [<Test>]
    member x.Mean() =
        let n = [0.;1.;2.;3.;4.;5.;6.]
        Assert.AreEqual(mean n, 3.)
        Assert.AreEqual(mean [], System.Double.NaN)
        Assert.AreEqual(mean [-1.;1.], 0.)
        Assert.AreEqual(mean [5.], 5.)
        
    [<Test>]
    member x.Median() =
        let n = [0.;1.;2.;3.;4.;5.;6.]
        Assert.AreEqual(median n, 3.)
        Assert.AreEqual(median [], System.Double.NaN)
        Assert.AreEqual(median [3.5], 3.5)
 
    [<Test>]
    member x.Sample() =
        let s = Seq.sample 5 a |> Seq.toList
        Assert.That(s |> Seq.forall (fun si -> Array.exists (fun ai -> ai = si) a))
        Assert.That( Seq.length s, Is.EqualTo(5) )
        
    [<Test>]
    member x.Memoize() =
        let n = ref 0
        let f a = incr n; a
        let memof = memoize f
        
        // The memoizer should only hit the underlying function once for each argument
        memof 2 |> ignore
        memof 3 |> ignore
        memof 2 |> ignore
        
        Assert.That(!n, Is.EqualTo(2))
        
        
    [<Test>]
    member x.Histogram() =
        let r = new Random(42)
        let n = 100
        let s = [ for i in {1..n} -> r.NextDouble() ]
        
        let (bins, binOcc) = histogram 10 s 
        // Check that all the points showed up in a bin
        Assert.That(n, Is.EqualTo(binOcc |> Seq.reduce (+)), "Missing some data points from the histogram")
        Assert.That(binOcc.Length, Is.EqualTo(10))


[<TestFixture>]
type TestMinimize() =

    [<Test>]
    member x.MakeSimplex() =
        let center_vec = Vector.ofArray [|0.;1.;2.;3.;4.;5.|]
        let length_vec = Vector.ofArray [|1.;2.;3.;4.;5.;6.|]
        
        let simplex = makeSimplex center_vec length_vec
        Assert.That(simplex.Length, Is.EqualTo((Vector.length center_vec) + 1))
        

    [<Test>]
    member x.NelderMead() =
        let rel_err (a:vector) (b:vector) = Vector.fold (fun a b -> Math.Max(a,b)) 0. ((a-b) .* (Vector.map (fun a -> 1./a) (a+b)))
        let c = Vector.create 3 -1.
        
        // A simple objective function (x-c)*x + 1
        let f x = Vector.sum ((x - c) .* x) + 1.
        let simplex = makeSimplex (Vector.create 3 -1.) (Vector.create 3 5.) 
        let ans = minimize f simplex (ObjFun(0.004)) id 100
                
        let diff = rel_err ans (Vector.create 3 1.5)
        
        Assert.That(diff, Is.LessThan(0.001))


[<TestFixture>]
type TestHmm() =
        
    // Evaluate some likelihoods and make sure they're reasonable
    [<Test>]
    member x.SimpleEval() =
        let t = "TTTTGCTAGTTGCTTGATCGCTCC"
        let r = "TTTTGTTTTGCTTCCTCCTCC"
        
        let model = makeModelFromEdnaParams (convDNA t) ednaParams.starter
        use eval = setupHmmEvaluators model (convDNA r)
        let likelihood = eval.p
        Assert.That(likelihood |> double |> Math.Log, Is.GreaterThan(-100.))
        Assert.That(likelihood |> double |> Math.Log, Is.LessThan(-2.))
        
        let t = "TTTTGCTAGTTGCTTGATCGCTCC"
        let r = "TTTTGTTTGACTTTGCTTCCTCCTCC"
        
        let model = makeModelFromEdnaParams (convDNA t) ednaParams.starter
        use eval = setupHmmEvaluators model (convDNA r)
        let likelihood = eval.p
        Assert.That(likelihood |> double |> Math.Log, Is.GreaterThan(-100.))
        Assert.That(likelihood |> double |> Math.Log, Is.LessThan(-2.))
        
    // Test that we can verify stochastic matrices
    [<Test>]
    member x.IsStochastic() =
        let good = Matrix.create 4 4 0.25
        Assert.That(isStochastic 0.001 good, Is.True, "4x4 w/ element 0.25")
        
        let bad = Matrix.create 4 4 0.2503
        Assert.That(isStochastic 0.001 bad, Is.False, "4x4 w/ element 0.2503")
        
        let x = 0.2
        let good2 = Matrix.init 4 4 (fun i j -> if i=j then 1. - x else x/3.)
        Assert.That(isStochastic 0.001 good2, Is.True)
        
    
    [<Test>]
    [<Category("Long")>]
    member x.AlignSweep() =
    
        let dumpCSV filename data =
            use f = new FileStream(filename, FileMode.Create)
            use w = new StreamWriter(f)
            
            w.WriteLine("enzyme, hmmfrac, acc, len, hmmtotal")
            
            // Write one line to the csv file
            let csvLine (name, frac, acc, len, hmmErr) =
                let sb = sprintf "%s, %f, %f, %f, %f" name frac acc len hmmErr
                w.WriteLine(sb)
                
            data |> Seq.iter csvLine
            
        let em q = Array.create 4 q
        let sp v (a:float[]) = Array.init (a.Length) (fun i -> v * a.[i])
             
                    
        let errs = [| 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.1; 1.2;|]
        let r = new Random(42)
        
        let lengths = [200; 400; 600; 800; 1000; 1200; 1400]
        let randLen() = lengths.[r.Next(lengths.Length)]
    
        
        let samp_acc errPars =
            let len = randLen()
        
            let tpl = [|1;2;3;4|] |> Seq.sampleWithRand rand len |> Seq.toArray
            
            let model = makeModelFromEdnaParams tpl errPars
            let read = sampleModel model
            
            let al = GlobalAlign.Align(iconvDNA tpl, iconvDNA read)
            
            (al.Accuracy, float(len))
            
            
        let alignN (pars, name) n fracErr =
            let p = pars fracErr
            
            Array.init n (fun _ -> 
                let (acc,len) = samp_acc (pars fracErr)
                (name, fracErr, acc, len, ednaParams.total_error p))
            
        
        let rows = [(samplePars, "453P"); ] |> // (parsbr, "Branchy")] |> 
                    Seq.collect (fun pars -> errs |> Seq.map (alignN pars 300)) |> Seq.concat        
        
        rows |> dumpCSV @"errByEnzyme4.csv"
        ()
    
    // Generate scaled error parameters
    member x.GenParams(e : float) =
        let sp v (a:float[]) = Array.init (a.Length) (fun i -> v * a.[i])
        
        let pars453 e : ednaParams =
            {   numbases = 4;
                dark = sp e [|0.05; 0.03; 0.04; 0.03|]
                insert = Array2D.init 4 4  (fun i j -> if i=j then 0.03 * e else 0.01 * e)
                merge = sp e [| 0.10; 0.1; 0.1; 0.1 |]
                miscall = Array2D.init 4 4 (fun i j ->
                    if i=j && i < 2 then 1. - 0.007 * e
                    else if i=j then 1. - 0.026 * e
                    else if i<2 && j<2 then 0.007 * e
                    else if i>=2 && j>=2 then 0.026 * e
                    else 0.); 
            }

        pars453 e
            

    // For trapping weirdness
    member x.LikelihoodPerBase() =
        let len = 1000

        let tpl = [|1;2;3;4|] |> Seq.sampleWithRand rand len |> Seq.toArray

        use ss = new FileStream(@"c:\users\pmarks\bad-edna.bin", FileMode.Create, FileAccess.Write, FileShare.None)
        let binf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter()

        let lpb e = 
            let model = makeModelFromEdnaParams tpl (x.GenParams e)
            let read() = sampleModel model       
                        
            let ab() =
                let r = read()
                use eval = setupHmmEvaluators model r
                let lpb = eval.p.Value / (float32 len)
            
                //let al = GlobalAlign.Align(iconvDNA tpl, iconvDNA r)

                if Single.IsNaN lpb || Math.Abs(lpb) > 10.0f then
                    printfn "%f" eval.p.Value
                    //printfn "%f" e
                    //printfn "%A" tpl
                    //printfn "%A" r 
                    binf.Serialize(ss, (tpl, r, model, eval.p.Value))

                    use eval2 = setupHmmEvaluators model r
                    printfn "%f" eval2.p.Value
                    printfn ""
                    let lpb2 = eval.p.Value / (float32 len)
                    binf.Serialize(ss, (tpl, r, model, eval2.p.Value))


                (1.0, lpb)

            
            let rr = Array.init 1000 (fun _ -> ab())
            (rr |> Array.map fst |> mean, rr |> Array.map (snd >> float) |> mean)

            
        let ee = [0.9; 0.95 ] |> List.map lpb
        printfn "%A" ee
        


    // Sample repeatedly from a small HMM to get empircal counts of the number of times we observe various output sequences.
    // Then compute the probability of each output seqeunce using the forward algo.  Compare the empirical to theoretical 
    // probabilities and make sure they jive.   
    [<Test>]
    member x.SampleCheck() =
    
        // Sample from the model and count the number of times each sequence is observed
        let sampleAndBin model n = 
            let m : Map<int[], int> ref = ref Map.empty
            
            
            for i = 1 to n do
                let s = sampleModel model
                m := match Map.tryFind s !m with
                        | Some(v) ->  Map.add s (v+1) !m
                        | None -> Map.add s 1 !m 
            done
            
            // Compute forward probabilty and empirical probability
            let checkSeq (key, value) =
                use eval = setupHmmEvaluators model key
                let p1 = eval.p |> double
                (p1, float(value) / float(n))

            (Map.toSeq !m) |> Seq.map checkSeq |> Seq.toList


        // Make a template sequence
        let intSeq = [|2;2;4;4;3;1|]   
        // Build the model
        let model = makeModelFromEdnaParams intSeq ednaParams.starter

        // Sample, and make sure it's close (binomial-wise) to the forward algorithm
        let n = 1e6
        let m = sampleAndBin model (int n)
        let ms = List.sortWith (fun (h1,_) (h2,_) -> if h1 > h2 then -1 else 1) m
        
        // Compute the sampling error and see how many sigmas of disagreement we have.
        let err_list = ms |> List.map (fun (prob,eprob) -> (prob, n * (eprob - prob) / Math.Sqrt(n * prob * (1.-prob)))) |> List.toArray        
        Assert.That(err_list |> Seq.forall (fun (prob, err) -> if prob > 1e4 && err > 4. then false else true), "Found HMM outlier")
    
 
    // Explore the impact of the Quiver scoreDiff parameter on Edna accuracy. 1e-8 seems to be the limit if you want exact results.
    [<Test>]
    member x.BandingSweep() =
        let accLevels = [| 1e-28; 1e-24; 1e-20; 1e-16; 1e-12; 1e-10; |]
        let t = convDNA "ACGTACGTAAGTCTACTACTACTACAACACAGAAGACATACGTCCGTAACACCGACCCAACCAAGATATACCCGGGTTACACGTGTGTACTACAAGATACGTCAGTACGAGATCACAGATACAGATATACCAACACGTGTACAGCTACGACACGACTACGACGACTACTACGACGGCGCGCATCGAGGAGAGGACCTCAGGCCTAGCGATCCACGAGCTCTACGACGACTTCAGACGCATACCGTACGACCGTACGGTACTTACGTACGTACGTACGTACAATGTACGGTACATACGACGTACGTTGCATGCAACGTCGTGCATACA";
        let fourCTrueModel = makeModelFromEdnaParams t ednaParams.starter 
        let dists = ednaParamsToTransEmDists (ednaParams.starter)        
        let r = sampleModel fourCTrueModel

        let doEval acc =
            
            let sw = new Stopwatch()
            
            sw.Start()
            let scoreDiff = float32 (Math.Abs(Math.Log(acc)))
            PacBio.Hmm.Recursions.scoreDiff := scoreDiff
            use eval = setupHmmEvaluators fourCTrueModel r
            let p = eval.p 
            sw.Stop()
                          
            (p, (sw.ElapsedMilliseconds, acc, scoreDiff))
            
            
        let res = accLevels |> Array.map doEval
        printfn "Template Length %d" t.Length
        printfn "%A" res
        
        let probs = res |> Array.map fst |> Array.map (fun v -> v.Value)
        printfn "log-probsL: %A" probs

        let errs = probs |> Seq.map(fun v -> float (Math.Abs(v-probs.[0]))) |> Seq.toArray        
        printfn "errors: %A" errs
        
        Assert.That( errs |> Seq.forall (fun v -> v < 1e-3), sprintf "Banding sweep found changes: %A" errs )          


    // Make sure alpha and beta agree
    [<Test>]
    member x.AlphaBetaCompare() =
        let t = "TGACAGACACATACATG"
        let model = makeModelFromEdnaParams (convDNA t) ednaParams.starter
        
        for i in 1 .. 100 do
            
            let r = sampleModel model
            
            printfn "read: %A:" r
            use eval = setupHmmEvaluators model r
            let alphaN = eval.p.Value
            let beta0 = (eval.beta).[0,0].Value
            
            printfn "alpha: %A, beta: %A" alphaN beta0            
            Assert.That(Math.Abs(alphaN - beta0) < 1e-2f)
        done        


// Run the EM algorithm and print out the sequence of likelihoods obtained.
// Verify that the likelihood increases monotonically.
let doEMParamterEstimation(reads, template, _guessParams, trueParams, trueModel) = 
    
    let guessParams = ednaParamsToTransEmDists _guessParams

    let getStatsUsingModel m = reads |> extractSufficientStatisticsSum m
    let reEstimateParams (p:transEmDists) = distsAndErrorsFromStats (getStatsUsingModel (makeModelFromTransEmDists 4 p template))
    
    let em m r = 
        use ev = (setupHmmEvaluators m r)
        ev.p.Value

    let readLikelihoods m = reads |> Seq.map (em m) |> Seq.reduce (+)

    let trueModel =makeModelFromTransEmDists 4 (ednaParamsToTransEmDists trueParams) template
    let trueLikelihood = readLikelihoods trueModel
    
    let updateModel (ntrys, p) = 
        if ntrys > 0 then 
            let newParams = reEstimateParams p
            let likelihood = readLikelihoods (makeModelFromTransEmDists 4 newParams.dists template)
            Some((newParams, likelihood), (ntrys-1, newParams.dists)) 
        else 
            None
    let emProgression = Seq.unfold updateModel (12, guessParams) |> Seq.toArray

    // Show what happened        
    Debug.WriteLine(sprintf "Likelihood of true model: %f" trueLikelihood)
    Debug.WriteLine("Likelihood progression:")

    let likelihoodString = sprintf "%A" (emProgression |> Array.map (fun (_,ll) -> ll))
    Debug.WriteLine(likelihoodString)
    
    // Make sure the EM increases the likelihood at every step, and that it gives a better model
    let (ordered, _) = 
        emProgression |> 
        Seq.fold (fun (ordered, last) (_, ll) -> if ll >= (last - 0.02f) && ordered then (true, ll) else (false, ll)) (true, Single.MinValue)

    let (bestP, bestLL) = emProgression |> Seq.reduce (fun (ap,al) (bp,bl) -> if al > bl then (ap,al) else (bp,bl))
    
        
    Debug.WriteLine(sprintf "True errors params: %A" (ednaParamsToTransEmDists trueParams))
    Debug.WriteLine(sprintf "Best model error params: %A" bestP.dists)
    Assert.That(ordered, Is.True,  "The EM did not converge monotonically: " + likelihoodString)
    Assert.That(bestLL >= trueLikelihood,
                (sprintf "The EM did not find a model at least as good as the true model: LL of found=%f, LL of true=%f" bestLL trueLikelihood))
    emProgression


let EMParamterEstimation(nreads, template, guessParams, trueParams) = 
    // The EM parameter estimation should generate a series of model parameters that
    // have monotonically increasing likelihood, and converge to a likelihood greater than
    // the true likelihood

    let trueModel = makeModelFromTransEmDists 4 (ednaParamsToTransEmDists trueParams) template

    let reads = [1..nreads] |> List.map (fun _ -> sampleModel trueModel)
    doEMParamterEstimation(reads, template, guessParams, trueParams, trueModel)

        
[<TestFixture>]
type ExpectationMaximization() =
    
    let t = convDNA "TGACAAGGCATCCTTCTACGAACTGGGTTAGTGGCCATGCAAGTGTGACCAACCCAGTGTGCGTAACAGACGTTTGGTTGTACACCCACAACGTGCATGTGCACCGCCGGCCGAAAAACTCTGTACGTACGTAGGTTTCCCACAGGTTAACCCGGTGTCACGTCTATAAACTT"


    let fourCTrueParams  = samplePars 1.0    
    let fourCGuessParams = samplePars 2.0        

    let fourCTrueModel = makeModelFromEdnaParams t fourCTrueParams
    
    let sm = [|0.000001; 0.000001; 0.000001; 0.000001|]
    
    let zeroErrorParams : ednaParams =
        { numbases=4;
            insert = Array2D.create 4 4 0.0001;
            miscall = Array2D.init 4 4 (fun i j ->  if i=j then 1. else 0. );         
            dark = sm;
            merge = sm;
         }
        
    let zeroErrorModel = makeModelFromEdnaParams t zeroErrorParams

    [<Test>]
    member x.ednaParamsToStatsAndBack() =
        let initParams = samplePars 0.8
        let dists = ednaParamsToTransEmDists initParams

        let pa = Array.create 4 1.0
        let da = Array.init 4 (fun i -> Array.create 5 1.0)
        let distErrs = { pStay = pa; pMerge = pa; stayDists = da; moveDists = da }
        let distsAndErr = { dists = dists; errors = distErrs }

        let rtParams = transEmDistsToEdnaParams distsAndErr
        let rtDist = ednaParamsToTransEmDists rtParams
        
        printfn "pars-start: %A" initParams
        printfn "pars-end: %A" rtParams


        printfn "dist-start: %A" dists
        printfn "pars-end: %A" rtDist


        let numbases = distsAndErr.dists.pStay.Length

        let distArray = transEmDists.toArray distsAndErr.dists
        let errorArray = transEmDists.toArray distsAndErr.errors

        if distArray |> Array.exists Double.IsNaN then 
            printfn "Got NaN in Param convert"

        let fmin (trialGenParamsVect : vector) = 
            let  trialGenParamsArray = trialGenParamsVect.InternalValues

            // lm may not respect the lb and ub that you give it, so fix the argument.
            let trialGenParamsArray = trialGenParamsArray |> Array.map (fun v -> Math.Max(MinP, Math.Min(MaxP, v)))
            let trialGenParams = ednaParams.ofArray numbases trialGenParamsArray

            let trialDistArray = transEmDists.toArray (ednaParamsToTransEmDists trialGenParams)

            let diffVect = Vector.init distArray.Length (fun i -> (distArray.[i] - trialDistArray.[i]) / errorArray.[i])
            diffVect


        let errSolution = fmin ((ednaParams.toArray rtParams) |> Vector.ofArray)

        let errTrue = fmin ((ednaParams.toArray initParams) |> Vector.ofArray)

        printfn "Solution norm: %A, Err: %A" errSolution (lm2.L2 errSolution)
        printfn "True norm: %A, Err: %A" errTrue (lm2.L2 errTrue)




    [<Test>]
    member x.ednaParamsRoundTrip() =
        // Make sure going back and forth between 4c error parameter and they flattened representation works
        let flat = ednaParams.toArray fourCTrueParams
        let rtp = (ednaParams.ofArray fourCTrueParams.numbases flat)
        Assert.AreEqual((ednaParams.toArray rtp), (ednaParams.toArray fourCTrueParams))   
       
    
    [<Test>]
    member x.transEmDistsRoundTrip() =
        let flat = transEmDists.toArray transEmDists.starter
        let roundTrip = transEmDists.ofArray flat
        Assert.AreEqual(flat, (transEmDists.toArray roundTrip))

    member x.DoTestFourCStats() =
        // Generate some reads, extract suff. stats from those reads,
        // and make sure that the 1C counts come back right
        let reads = [1..500] |> List.map (fun _ -> sampleModel fourCTrueModel)
        //let tempStats = reads |> Seq.pmap (extractSufficientStatistics fourCTrueModel) |> Seq.toArray
        //let empiricalStats = tempStats |> Seq.reduce (+)
        
        let suffStats = reads |> extractSufficientStatisticsSum fourCTrueModel

        let distsAndErrs = distsAndErrorsFromStats suffStats
        let observedDists = transEmDists.toArray distsAndErrs.dists
        let observedErrs = transEmDists.toArray distsAndErrs.errors

        let returnParams = transEmDistsToEdnaParams distsAndErrs
        let expectedDists = transEmDists.toArray (ednaParamsToTransEmDists fourCTrueParams)
        
        let nbases = reads.Length * t.Length
        // Determining what is 'significantly different' here is a bit tough -- basically all the numbers are 
        // probilities, and the number of examples you saw is the template length / 4 * the probability
        let sigma = 
            Seq.zip3 expectedDists observedDists observedErrs |> 
            Seq.map (fun (expected,obs, obsErr) -> 
                if expected < 1e-4 then 0.0
                else
                    (obs-expected) / obsErr) |> 
            Seq.toArray
            
        printfn "Sufficient Statistics"
        printfn "%A" suffStats

        printfn "True Parameters"
        printfn "%A" fourCTrueParams
        
        printfn "Fitted Parameters"
        printfn "%A" returnParams

        printfn "Observed Probabilites"
        printfn "%A" observedDists
        printfn "Expected Probabilities"
        printfn "%A" expectedDists
        printfn "Got sigmas of suff. stats difference:\n %A" sigma
        
        //(sigma |> Seq.filter Double.IsNaN |> Seq.length < 10) //&&
        (sigma |> Seq.filter (Double.IsNaN >> not) |> Seq.forall (fun e -> Math.Abs(e)<6.))
        
        
    [<Test>]
    member x.fourCSufficientStatistics() =        
        let rec tryTest n = 
            if n = 0 then Assert.Fail "Error counts were significantly different from nominal"
            if x.DoTestFourCStats() then 
                ()
            else
                tryTest (n-1)
                
        tryTest 3

        

    [<Test>]
    member x.fourCZeroErrorSufficientStatistics() =
        // Generate some reads, extract suff. stats from those reads,
        // and make sure that the 1C counts come back right
        let nreads = 1000
        let reads = [1..nreads] |> List.map (fun _ -> sampleModel zeroErrorModel)
        let empiricalStats = reads |> extractSufficientStatisticsSum zeroErrorModel

        let distsAndErrs = distsAndErrorsFromStats empiricalStats
        let counts = transEmDists.toArray distsAndErrs.dists

        let returnParams = transEmDistsToEdnaParams distsAndErrs

        //let (returnParams, counts) = fourCEstimateParams empiricalStats fourCGuessParams
        let expectedCounts = transEmDists.toArray distsAndErrs.dists
        
        // Determining what is 'significantly different' here is a bit tough -- basically all the numbers are 
        // probilities, and the number of examples you saw is the template length / 4 * the probability
        let nbases = nreads * t.Length
      
        printfn "True Params"
        printfn "%A" zeroErrorParams
        
        printfn "New Params"
        printfn "%A" returnParams
        
        printfn "Empirical Counts"
        printfn "%A" empiricalStats
        printfn "Empirical Probs"
        printfn "%A" counts
        printfn "Expected Probs"
        printfn "%A" expectedCounts
        
        Assert.That(returnParams.merge |> Seq.forall (fun e -> e < 0.01));
        Assert.That(returnParams.dark  |> Seq.forall (fun e -> e < 0.01));     
        {0..3} |> Seq.map (fun i -> returnParams.miscall.[i, i]) |> Seq.forall (fun e -> e > 0.995) |> Assert.That
        

    [<Test>]
    member x.fourCEMParameterEstimationBig() = 
        EMParamterEstimation(300, t, fourCGuessParams, fourCTrueParams) |>
        ignore
        
    [<Test>]
    member x.fourCEMParameterEstimationSmall() =         
        EMParamterEstimation(30, t, fourCGuessParams, fourCTrueParams) |>
        ignore
        
    
    [<Test>]
    member x.fourCEMParameterEstimationTrueStart() = 

        let nreads = 100
        let reads = [1..nreads] |> List.map (fun _ -> sampleModel fourCTrueModel)
        
        let emProg = doEMParamterEstimation(reads, t, fourCTrueParams, fourCTrueParams, fourCTrueModel)
        ()

    [<Test>]
    member x.fourCZeroErrors() =     
    
        let nreads = 100
        let reads = [1..nreads] |> List.map (fun _ -> sampleModel zeroErrorModel)

        let _ = doEMParamterEstimation(reads, t, fourCGuessParams, zeroErrorParams, zeroErrorModel)
        ()

    // Run EM on the same data from a bunch of different starting points and make sure we get the same answer
    [<Test>]
    member x.fourCEMParameterEstimationDifferentStarts() =
        let g1 : ednaParams =  
            { numbases=4; 
            merge=[|0.2; 0.1; 0.1; 0.05|];
            miscall = Array2D.init 4 4 (fun i j ->
                 if i=j then 1. - 0.03 
                 else if i<2 && j<2 then 0.03
                 else if i>=2 && j>=2 then 0.03
                 else 0.);
            insert = Array2D.init 4 4 (fun i j -> if i=j then 0.05  else 0.03) ;
            dark = [|0.05; 0.05; 0.05; 0.05|] }          

        let g2 : ednaParams =  
            { numbases=4; 
            merge=[|0.1; 0.3; 0.05; 0.2|];
            miscall = Array2D.init 4 4 (fun i j ->
                 if i=j then 1. - 0.01 
                 else if i<2 && j<2 then 0.01
                 else if i>=2 && j>=2 then 0.01
                 else 0.);
            insert = Array2D.init 4 4 (fun i j -> if i=j then 0.02  else 0.03) ;
            dark = [|0.15; 0.05; 0.15; 0.05|] }          

        let g3 : ednaParams =  
            { numbases=4; 
            merge=[|0.2; 0.1; 0.1; 0.05|];
            miscall = Array2D.init 4 4 (fun i j ->
                 if i=j then 1. - 0.03 
                 else if i<2 && j<2 then 0.03
                 else if i>=2 && j>=2 then 0.03
                 else 0.);
            insert = Array2D.init 4 4 (fun i j -> if i=j then 0.05  else 0.03) ;
            dark = [|0.05; 0.05; 0.05; 0.05|] }          

        let g4 : ednaParams =  
            { numbases=4; 
            merge=[|0.03; 0.01; 0.01; 0.05|];
            miscall = Array2D.init 4 4 (fun i j ->
                 if i=j then 1. - 0.01 
                 else if i<2 && j<2 then 0.01
                 else if i>=2 && j>=2 then 0.01
                 else 0.);
            insert = Array2D.init 4 4 (fun i j -> if i=j then 0.02  else 0.03) ;
            dark = [|0.08; 0.02; 0.05; 0.07|] }          

                          
        let nreads = 150
        let reads = [1..nreads] |> List.map (fun _ -> sampleModel fourCTrueModel)

        let r1 = doEMParamterEstimation(reads, t, g1, fourCTrueParams, fourCTrueModel)
        let r2 = doEMParamterEstimation(reads, t, g2, fourCTrueParams, fourCTrueModel)
        let r3 = doEMParamterEstimation(reads, t, g3, fourCTrueParams, fourCTrueModel)
        let r4 = doEMParamterEstimation(reads, t, g4, fourCTrueParams, fourCTrueModel)

        let llr b = b |> Array.map snd |> Array.sum

        let llFinal = [| llr r1; llr r2; llr r3; llr r4 |]

        printf "Final likelihoods for different starting points: %A" llFinal 
        ()

    [<Test>]
    member x.EMBandingSweep() =
        let accLevels = [| 1e-8; 1e-6 |] |> Array.map (fun v -> new LFloat(v))
        let t = convDNA "ACGTACGGTAAGTCTTACTACTACTACAACACAGAAGGACATACGTCCGTTAAGCACCGACCCAACCAAGATATACCCGGGTTACACGTGTGTACTTACAAGATA";
        let fourCTrueModel = makeModelFromEdnaParams t ednaParams.starter        
        let r = sampleModel fourCTrueModel        
        
        let doEval acc =
        
            let sw = new Stopwatch()            
            sw.Start()
            let ss = extractSufficientStatistics fourCTrueModel r
            sw.Stop()

            (ss, sw.ElapsedMilliseconds)            
            
        let res = accLevels |> Array.map doEval
        printfn "Template Length %d" t.Length
        printfn "%A" res

[<TestFixture>]
type LevenbergMarquardt() =

    // Try to find the parameters of a simple function from some noisy sampling of it.
    [<Test>]
    member x.TrigMinimize() =
    
        let r = new Random(42)
        
        // The parameters to generate the fitting data
        let trueParams = Vector.ofArray [| 10.; 20.; 6.; 50. |]
        // The starting fit parameters
        let guessParams = Vector.ofArray [| 8.; 18.; 4.; 44.|]
        
        // The function taking the parameters and an x value and giving a y value
        let f (p : Vector<float>) x = p.[0] * Math.Cos(x * p.[1]) + p.[2] * x + p.[3];
        
        // The x values
        let xvals = [|0..10000|] |> Array.map(fun i -> float(i)/10000.)
        
        // The true y values
        let y = xvals |> Array.map (fun x -> (f trueParams x)) //)
        
        let residuals pars = xvals |> Array.mapi (fun i x -> y.[i] - (f pars x)) |> Vector.ofArray
        
        let ansParams = lm2.lm residuals guessParams (lm2.OptimizationOptions.standard 100)        
        let errors = Seq.zip (Vector.toArray trueParams) (Vector.toArray ansParams.p) |> Seq.map (fun (t,a) -> (t-a) / t)
        
        printfn "True Parameters:\n%A" trueParams
        printfn "Got solution:\n%A" ansParams
        Assert.That(errors |> Seq.forall (fun e -> e < 1e-3), "LevenbergMarquardt didn't converge on solution");     
        

    // Try to find the parameters of a simple function from some noisy sampling of it.
    [<Test>]
    member x.PolyMinimize() =
        
        // The parameters to generate the fitting data
        let trueParams = Vector.ofArray [| 2.; 3.; 1.; 0.1 |]
        // The starting fit parameters
        let guessParams = Vector.ofArray [| 1.; 1.; 1.; 1.|]
        
        // The function taking the parameters and an x value and giving a y value
        let f (p : vector) x = p.[0] + p.[1] * x + p.[2] * x**2. + p.[3] * x**3.
        
        // The x values
        let xvals = [|-20 .. 20|] |> Array.map(fun i -> float(i)/5.)
        
        // The true y values
        let y = xvals |> Array.map (f trueParams)
        
        let residuals pars = xvals |> Array.mapi (fun i x -> y.[i] - (f pars x)) |> Vector.ofArray
        
        
        let ansParams = lm2.lm residuals guessParams (lm2.OptimizationOptions.standard 100)        
        let errors = Seq.zip (Vector.toArray trueParams) (Vector.toArray ansParams.p) |> Seq.map (fun (t,a) -> (t-a) / t)
        
        printfn "True Parameters:\n%A" trueParams
        printfn "Got solution:\n%A" ansParams
        Assert.That(errors |> Seq.forall (fun e -> e < 1e-3), "LevenbergMarquardt didn't converge on solution");      
        
    // Make sure our box/linear constrained minimizer works as it should 
    [<Test>]
    member x.ConstrainedMinimize() =        
                
        // The parameters to generate the fitting data
        let trueParams = Vector.ofArray [| 2.; 3.; 1.; 0.1 |]
        // The starting fit parameters
        let guessParams = Vector.ofArray [| 6.1; 0.; 0.; 0.|]
        
        // The function taking the parameters and an x value and giving a y value
        let f p x = 
            Assert.That( Math.Abs((p |> Vector.toArray |> Seq.reduce (+)) - 6.1) <  0.0001 )
            p.[0] + p.[1] * x + p.[2] * x**2. +  p.[3] * x**3.
        
        // The x values
        let xvals = [|-20 .. 20|] |> Array.map(fun i -> float(i)/5.)
        
        // The true y values
        let y = xvals |> Array.map (f trueParams)
        
        let residuals pars = xvals |> Array.mapi (fun i x -> y.[i] - (f pars x)) |> Vector.ofArray
        
        let lb = Vector.create 4 -1.0 
        let ub = Vector.create 4 7.0
        let A = Array2D.create 1 4 1.0 |> Matrix.ofArray2D
        let b = [| 6.1 |] |> Vector.ofArray
        
        let ansParams = lm2.lmBoxConstrains (Some lb) (Some ub) A b residuals guessParams (lm2.OptimizationOptions.standard 100)
        
        let errors = Seq.zip (Vector.toArray trueParams) (Vector.toArray ansParams.p) |> Seq.map (fun (t,a) -> (t-a) / t)
        
        printfn "True Parameters:\n%A" trueParams
        printfn "Got solution:\n%A" ansParams
        Assert.That(errors |> Seq.forall (fun e -> e < 1e-3), "LevenbergMarquardt didn't converge on solution");  

