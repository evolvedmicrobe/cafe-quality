#light
module PacBio.Hmm.Minimize

open System
open System.Collections.Generic
open Microsoft.FSharp.Math
open System.Diagnostics


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



// Easy printing routine
let printl x = Console.WriteLine (sprintf "%A" x)

// Stupid simplex generator -- take a guess an go out along each vector u dot scale, to make n points and use center as n+1th
let makeSimplex (center : vector) (scale : vector) =
    let uv i = Vector.init center.Length (fun j -> if j=i then 1. else 0.)
    let mkPt i = if i=0 then center else center + ((uv (i-1)) .* scale)
    Array.init (center.Length + 1) mkPt 
    

let gt a1 (a2 : vector) = Vector.mapi (fun x v1 -> if v1 > a2.[x] then 1. else 0.) a1
let vcom f (v1:vector) (v2:vector) = Vector.mapi (fun i e1 -> f(e1, v2.[i])) v1
let vmax = vcom Math.Max
let vmin = vcom Math.Min

type tol = ObjFun of float | Argument of vector

let minimize (func : vector -> float) (initsimplex : vector array) tol constrain maxFuncIters =
    // Nedler-Mead simplex minimization. func is the objective function to minimize, it must take a n-dimensional vector and returns a float. initsimplex is the 
    // initial simplex to use -- it is a array of n+1 n-dimension vectors, which hopefully are in the region, or bound the minimum you are looking for. tol is the 
    // relative error in the objective function you can tolerate.
    // Returns a n-dimensional vector containing the best point found.
    
    let (f, funcIters) = nmemoize func
    let ndims = initsimplex.[0].NumRows
    let psum simplex = Vector.init ndims (fun i -> simplex |> Seq.fold (fun c (pt : vector) -> c + pt.[i]) 0.) 
    
    let findBounds (simplex : vector array) =
        let idx_arr = [|0 .. simplex.Length - 1|]
        
        let sortf i1 i2 =
            if (f simplex.[i1]) > (f simplex.[i2]) then
                -1 
            else if (f simplex.[i1]) = (f simplex.[i2]) then 
                0 
            else 
                1
        
        let idx_arr = Array.sortWith sortf [|0 .. simplex.Length - 1|]
        (idx_arr.[0], idx_arr.[1], idx_arr.[idx_arr.Length - 1])            
            
    let amotry (simplex : vector[]) ihi fac =
        let fac1 = (1.0 - fac) / float(ndims)
        let fac2 = fac1 - fac
        let ptry = constrain ((psum simplex) * fac1 - simplex.[ihi] * fac2)
        
        match f ptry with
            | x when x < f simplex.[ihi] -> simplex.[ihi] <-  ptry ; x
            | x -> x
    
    let contract simplex = 
        let (_,_, ilo) = findBounds simplex
        let shrink pt = (Vector.create ndims 0.5) .* (pt + simplex.[ilo])
        Array.iteri (fun i pt -> if i = ilo then () else simplex.[i] <- constrain (shrink simplex.[i])) simplex
        
    
    let rec amoeba simplex =
        let (yhi, ynhi, ylo) = findBounds simplex
        let update simplex = 
            match (amotry simplex yhi -1. , f simplex.[yhi]) with
                | x, _ when x <= f simplex.[ylo] -> ignore(amotry simplex yhi 2.)
                | x, ysave when x >= f simplex.[ynhi] -> if (amotry simplex yhi 0.5) > ysave then contract simplex else ()
                | _ -> ()
        
        let ftol = 2.0*abs((f simplex.[yhi]) - (f simplex.[ylo]))/(abs ((f simplex.[yhi]) + (f simplex.[ylo]) + 1e-10) )
        let avgpt = (1./float(simplex.Length)) * (simplex |> Seq.reduce (+)) 
        let stol = simplex |> Seq.map (fun s -> s - avgpt) |> Seq.map (Vector.map abs) |> Seq.reduce vmax
        let r =  (f simplex.[yhi])
        
        let tolmet = match tol with
                        | ObjFun(t) -> t > ftol
                        | Argument(t) -> Vector.forall (fun x -> x > 0.) (gt t stol)
        
        match (tolmet, funcIters()) with
            | (false, n) when n < maxFuncIters -> update simplex; amoeba simplex
            | (false, _) -> Debug.WriteLine("minimization terminated because MaxIterations was exceeded"); simplex.[ylo]
            | _ -> simplex.[ylo]
            
    amoeba initsimplex