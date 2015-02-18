#light
module PacBio.Hmm.Utils

open System
open System.Collections.Generic
open System.Reflection
open Microsoft.FSharp.Math
open System.Linq
open System.Threading
open System.Threading.Tasks


// Identity function
let id a = a

// Set up some constraints.
let vcom f (v1:vector) (v2:vector) = Vector.mapi (fun i e1 -> f(e1, v2.[i])) v1
let vmax = vcom Math.Max
let vmin = vcom Math.Min
let vs = Vector.ofArray [| 0.001; 0.001; 0.001; 0.001 |]
let vl = Vector.ofArray [| 0.999; 0.999; 0.999; 0.999 |]
let pconstrain v = 
    let v = vmin vl (vmax v vs)
    v

let memoize f = 
    // Simple memoizer
    let cache = Dictionary<_, _>()
    let ev_func = fun x -> 
        let ok,res = cache.TryGetValue(x)
        if ok then res 
        else let res = f x
             cache.[x] <- res
             res
    ev_func


let nmemoize f = 
    // Simple memoizer with 
    // function invocation counting
    let n = ref 0
    let cache = Dictionary<_, _>()
    let ev_func = fun x -> 
        let ok,res = cache.TryGetValue(x)
        if ok then res 
        else let res = f x
             incr n
             cache.[x] <- res
             res
    (ev_func, (fun () -> !n))

module Seq =
    let r = new Random(42)
    // Sample with replacement
    // Be careful to cache the result -- otherwise people will get a different result every time the enumerate the returned seq
    let sample n s = 
        let a = s |> Seq.toArray
        let len = a.Length
        if len = 0 then
            Seq.empty
        else
            let l = [ for i in {1..n} -> a.[r.Next(len)] ]
            l |> Seq.ofList
        
    let sampleWithRand (rand : Random) n s = 
        let a = s |> Seq.toArray
        let len = a.Length
        if len = 0 then
            Seq.empty
        else
            let l = [ for i in {1..n} -> a.[rand.Next(len)] ]
            l |> Seq.ofList

    let pop s = Seq.find (fun _ -> true) s
    
    let chooseWhile f =
        Seq.initInfinite (fun _ -> f()) |> Seq.takeWhile Option.isSome |> Seq.choose id


module List =
    // Skip the first n items of the list
    let rec skip n l =
        match n with
            | n when n > 0 -> skip (n-1) (List.tail l)
            | _ -> l

module Array =
    let rec Slice (s : #seq<int>) (a : 'a[])  =
        let na = Array.zeroCreate (Seq.length s)
        s |> Seq.iteri (fun si v -> na.[si] <- a.[v])
        na

let stats s =
    match Seq.length s with
        | 0 -> (Double.NaN, Double.NaN)
        | _ -> 
            let (sum, sq_sum, n) = s |> Seq.fold (fun (sum, sq_sum, n) b -> (sum+b, sq_sum+b*b, n+1)) (0., 0., 0)
            (sum/float(n), Math.Sqrt(sq_sum-sum*sum)/float(n))
            

let mean s = 
    match Seq.length s with
        | 0 -> System.Double.NaN
        | _ -> (Seq.fold (+) 0. s) / float(Seq.length s)

// Find the median of a seq of floats
let median (s : #seq<float>) = PacBio.Utils.Stats.QuickMedian(Seq.toArray s)

let histogram nbins (data : #seq<float>) =
    let arr = data |> Seq.toArray
    
    if arr.Length = 0 then
        (Array.create 0 0., Array.create 0 0)
    else
    
        let min = 0.1 //Stats.Percentile(arr, 0.1);
        let max = 0.9 //Stats.Percentile(arr, 80.);
        
        let h = Math.Max(0.01, (max - min) / float(nbins))
        let binBots = {0 .. nbins-1} |> Seq.map (fun v -> h * float(v) + min) |> Array.ofSeq
        let bins = Array.create nbins 0
        
        let putInBin (v : float) = 
            let b = Math.Max(0, Math.Min(int(Math.Floor((v - min) / h)), nbins-1))
            bins.[b] <- bins.[b] + 1
        
        data |> Seq.iter putInBin
        
        (binBots, bins)
    
    
let r = new System.Random()

// Easy printing routine
let printl x = Console.WriteLine(sprintf "%A" x) 
// Membership test 
let ismember set item  = Seq.exists (fun e -> e = item) set    

let randsample rng n = 
    let rec samp rng = 
        match r.Next(rng) with
            | 0 -> samp rng
            | v -> v
    Array.init n (fun _ -> samp rng)


let LogEpsilon = Math.Log(System.Double.Epsilon)
let LogEpsilon2 = LogEpsilon / 2.
let ExpansionMin = LogEpsilon / 4.

let MinP = Math.Sqrt(Double.Epsilon)
let MaxP = 1.0 - 0.0000001

// Add an extension method for doing underflow checked log-space addition
type System.Math with
    static member inline LogAdd(a : float, b : float) =
        let (mx,mn) = if a > b then (a,b) else (b,a) 
        
        match (mn-mx) with
            | d when d < LogEpsilon -> mx
            | d when d < LogEpsilon2 -> mx + Math.Exp(d)
            | d when d < ExpansionMin -> 
                let r = Math.Exp(d)
                let rsq = r*r
                mx + r - rsq/2. + rsq*r/3.
            | d when System.Double.IsNaN(d) = false -> mx + Math.Log(1. + Math.Exp(d))
            | _ -> mx


let convDNA (s : string) = s.ToCharArray() |> Array.map (function 'A' -> 1 | 'C' -> 2 | 'G' -> 3 | 'T' -> 4 | _ -> 0)
let convT = convDNA

let iconvDNA (i : int[]) = i |> Array.map (function 1 -> 'A' | 2 -> 'C' | 3 -> 'G' | 4 -> 'T' | _ -> 'N') |> (fun c -> new String(c))
