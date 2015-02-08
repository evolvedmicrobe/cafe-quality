module PacBio.Analysis.Models
(*
open System
open PacBio.Common.Numeric
open PacBio.Analysis.DataWrap
open PacBio.Analysis.Data
open Distributions

// Define a more realistic model that takes into account the frame rate of the camera and
// the fact that you can only acquire an integer number of frames.

let sample_missing tau_f photons_f baseline (thresholds : Map<int,double>) cameraGain (camPdf : Distributions.CameraPDFs) (r:Random) =   
    let rexp(mu) = -mu * Math.Log(r.NextDouble())
    
    let tstart = r.NextDouble()
    let pw = rexp(tau_f)
    let tend = tstart + pw
    let tend_long = Math.Ceiling(tend) 
    let tend_shrt = Math.Floor(tend)
    
    let longestTimeScale = thresholds |> Seq.map (fun p -> p.Key) |> Seq.reduce max    
    let bgFrames v = Math.Max(1, Math.Min(longestTimeScale, int(v)))
    let detectedPulseTime t = Math.Max(0., Math.Min(float(longestTimeScale), t))
    
    // Signal photons
    let s1 = detectedPulseTime pw
    let s2 = detectedPulseTime (tend_shrt - tstart)
    let s3 = detectedPulseTime (tend      - 1.)
    let s4 = detectedPulseTime (tend_shrt - 1.)   
    
    // Background photons
    let b1 = bgFrames tend_long
    let b2 = bgFrames (Math.Max(1., tend_shrt))
    let b3 = bgFrames (Math.Max(1., tend_long-1.))
    let b4 = bgFrames (Math.Max(1., tend_shrt-1.))
    
    // Sigma correction - convert a threshold from the signal distribution to the
    // background distribution   
    let sigc s n = Math.Sqrt(2. * n) / Math.Sqrt(2. * (s + n))
    // The signal-to-noise of the pulse, including the signal variance
    let snr s n = s / Math.Sqrt(2. * (s + n))    
    // compute the missing probability
    let snr(s,nbg) = (s * photons_f + float(nbg) * baseline) / thresholds.[nbg]
    // Which slicing gives us the best chance of detecting the pulse?
    let sbgPairs = [|(s1,b1); (s2,b2); (s3,b3); (s4,b4)|]
    let mps = sbgPairs |> Array.map snr 
    let (s,b) = sbgPairs.[Stats.IMax mps]
    let threshold = thresholds.[b]
    let signal = s * photons_f + float(b) * baseline
    
    let nsample = 12
    let ndetect = {1..nsample} |> Seq.filter (fun _ -> camPdf.EMCCDSample(signal, cameraGain) > threshold) |> Seq.length
    
    1.0 - float(ndetect) / float(nsample)
    
    
let meanMissing n sigma tau_s photons_f baseline_f frame_time emccdgain =    
    if Double.IsNaN(photons_f) then Double.NaN else
        // Get the EMCCD corrected threshold in counts for the given baseline and sigma
        let camPdf = new Distributions.CameraPDFs();
        
        let r = new Random()
        
        let getThresh bg = camPdf.EMCCDThreshold(sigma, bg, emccdgain)
        let thresholds = {1..4} |> Seq.map (fun i -> (i,getThresh (baseline_f * float(i)))) |> Map.ofSeq
        Array.init n (fun _ -> sample_missing (tau_s/frame_time) photons_f baseline_f thresholds emccdgain camPdf r) |> Seq.average
    
    
let brightness (p : Pulse[]) = 
    [| for i in 0..3  do yield p |> Seq.filter (fun p -> p.Channel = i) |> Seq.map (fun p -> p.PkMid) |> Stats.NaNMedian |]

    
let analyticalMissing (sigmas : float[]) (tau_s : float[]) (tr : Trace) (p : Pulse[]) =
    // Pull out the aligning pulses
    let br = brightness p
    
    // Use the semi-analytical model of missing fraction to compute the expected missing rate
    // for the brightness and BG of the trace
    let missing ch =
        try
            meanMissing 1000 sigmas.[ch] tau_s.[ch] (br.[ch]) (float tr.BaselineBias.[ch]) (1./tr.Scan.FrameRate) (float tr.Scan.Metadata.CameraGain)
        with
            _ -> Double.NaN
    
    let r = [| for i in 0..3 do yield missing i |]
    if r |> Seq.exists (fun v -> v < 0.01) && System.Diagnostics.Debugger.IsAttached then
        System.Diagnostics.Debugger.Break()
        
    r
    
// Compute the stick rate we should observe from stochastic false positives at a given sigma level    
let stochasticFalsePositive sigma pulseRate fps meanDt meanIpd =
    let pFrame = normaldistr.normaldistribution(-sigma)
    let sticksPerPulse = fps * pFrame / pulseRate * (meanIpd / (meanDt + meanIpd))
    
    sticksPerPulse
    
let analyticalFalsePositive (sigmas : float[]) (tau_s : float[]) (ipd_s : float[]) fps =

    // Use the semi-analytical model of missing fraction to compute the expected missing rate
    // for the brightness and BG of the trace
    let pulseRate = 1.0 / ((ipd_s |> Seq.average) + (tau_s |> Seq.average))
    
    let fp ch =
        stochasticFalsePositive sigmas.[ch] pulseRate fps (Seq.average tau_s) (Seq.average ipd_s)
    
    // There's only a 1/4 chance the next base will be the same as the current base.    
    // We are neglecting the possibility of different colour branching.
    [| for i in 0..3 do yield (fp i) |]    



let sampleMerging ipdTauFrames pkMid baseline downThreshold fps cameraGain (camPdf : Distributions.CameraPDFs) (r : Random) =
    
    let rexp(mu) = -mu * Math.Log(r.NextDouble())
    
    let nf = rexp(ipdTauFrames)

    if nf > 1.0 then 
        // If it's off for more than one frame then it's not a merge        
        0.0 
    else
        // If it's off for less than one frame, then compute the chances of the intensity
        // failing to drop below the downThreshold.
        // *** Currently this analysis assumes that the IPD all falls inside one frame.  A more careful analysis 
        // *** that considers IPD crossing frame boundaries is required.
        
        // Number of photons in this frame
        let signal = baseline + pkMid * (1.0 - nf)
        
        let nsample = 20    
        // Compute the probability that this exceeds the down threshold
        let ndetect = {1..nsample} |> Seq.filter (fun _ -> camPdf.EMCCDSample(signal, cameraGain) > downThreshold) |> Seq.length    
        float(ndetect) / float(nsample)   
   
    
let meanMerging n downSigma medianPkmid baseline ipdTau fps emccdgain =    
    let camPdf = new Distributions.CameraPDFs();    
    let getThresh bg = camPdf.EMCCDThreshold(downSigma, baseline, emccdgain)
    let r = new Random()
    let downThreshold = getThresh baseline
        
    Array.init n (fun _ -> sampleMerging (ipdTau * fps) medianPkmid baseline downThreshold fps emccdgain camPdf r) |> Seq.average
    
    
    
let analyticalMerging (downSigmas : float[]) (ipd_s : float[]) (tr : Trace) (p : Pulse[]) =
    // Pull out the aligning pulses
    let br = brightness p
    
    // Use the semi-analytical model of missing fraction to compute the expected missing rate
    // for the brightness and BG of the trace
    let merging ch =
        try
            meanMerging 1000 downSigmas.[ch] (br.[ch]) (float tr.BaselineBias.[ch]) ipd_s.[ch] tr.Scan.FrameRate (float tr.Scan.Metadata.CameraGain)
        with
            _ -> Double.NaN
    
    // There's only a 1/4 chance the next base will be the same as the current base.    
    // We are neglecting the possibility of different colour branching.
    [| for i in 0..3 do yield (merging i) / 4.0 |]

*)