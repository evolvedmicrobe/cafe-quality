module lm2

open MathNet.Numerics.LinearAlgebra.Factorization

open MathNet.Numerics.LinearAlgebra.Double
open MathNet.Numerics.LinearAlgebra
open System

let L2 = Math.Vector.fold (fun acc v -> acc + v**2.) 0.0
type OptimizationOptions = {
    Delta : float;
    Eps1 : float;
    Eps2 : float;
    Eps3 : float;
    Tau : float;
    MaxIters : int;
    }

type OptimizationOptions with
    static member standard(maxIters) = {
        Tau = 1e-3;
        Eps1 = 1e-17;
        Eps2 = 1e-17;
        Eps3 = 1e-17;
        Delta = 1e-6;
        MaxIters = maxIters
    }

MathNet.Numerics.Control.DisableParallelization <- true

// Find the Y and Z such that x = Y * b + Z * x_z  and Ax=b for all x_z
let elim (A : matrix) (b : vector) =
    let m = A.NumRows // number of constraints
    let n = A.NumCols // number of variables

    // Find 'rank-revealing' QR decomposition of A^T
    let dm = DenseMatrix.OfArray(A.Transpose.ToArray2D())
    MathNet.Numerics.LinearAlgebra.Factorization.Dens
    let qr = new DenseQR(dm)
    let _ = qr.Solve(dm) 

    let r = qr.R
    let eps = 1e-10
    let rank = m

    let rInvTrans = r.SubMatrix(0, rank, 0, rank).Inverse().Transpose()

    if not qr.IsFullRank then 
        failwith("Constraint matrix must be full rank")

    let q = qr.Q
    let y = q.SubMatrix(0, n, 0, rank) * rInvTrans
    let Y = y.ToArray() |> Matrix.ofArray2D

    let c = Y * b
    let Z = q.SubMatrix(0, n, rank, n-rank).ToArray() |> Matrix.ofArray2D
    (Y, Z, c)



let jacobian (f : vector -> vector) (p : vector) (hx : vector) delta =
    //let jac = Matrix.create (hx.Length) (p.Length) 0.0
    let jac = new MathNet.Numerics.LinearAlgebra.Double.DenseMatrix(hx.Length, p.Length)
    
    for j in 0..(p.Length - 1) do
        let d = max delta (p.[j] * 1e-4)
        let tmp = p.[j]

        p.[j] <- p.[j] + d
        let hxx = f p
        p.[j] <- tmp

        for i in 0..(hx.Length - 1) do
            jac.[i,j] <- (hxx.[i] - hx.[i]) / d

    jac

(*
// WIP in Schnabel's backtracking line search method
let lineSearch (f : vector -> vector) (x : vector) (p : vector) (grad : vector) =

    let lambda = 1.

    let xpls = x + lambda * p

    let hx = f xpls

    // Loop
    let loop lambda k =

        let hx_L2 = L2 hx
        let fpls = 0.5 * hx_L2

        let kNew =
            let tLambda =  
                if k = 0 then 
                    -lambda * slp / (fpls - f - slp) * 0.2 
                else
                    let t1 = fpls - f - lambda * slp
                    let t2 = pfpls - f - plmbda * slp
                    let t3 = 1. / (lambda - plmabda);
                    let a3 = 3. * t3 * (t1 / (lambda**2) - t2 / (plambda**2))
                    let b = t3 * (t2 * lambda / plambda**2 - t1 * plambda / lambda**2)
                    let disc = b**2 - a3 * slp

                    let tl = if (disc > b * b) then
                        //  only one positive critical point, must be minimum */
                        (-b + if a3 < 0 then -Math.Sqrt(disc) else Math.Sqrt(disc)) / a3;
                    else
                        ///* both critical points positive, first is minimum */
                        (-b + if a3 < 0 then Math.Sqrt(disc) else -Math.Sqrt(disc)) / a3;

                    let tlambda = min(tl, lambda * 0.5)



    // Check for solution
    //if fpls <= f0 + slp * alpha * lambda then
*)

type constraints = { ub : vector option; lb : vector option }
type constraints with
    member x.Clip (v : vector) = 
        let clip i (v : float) =
            let c1 = match x.ub with None -> v | Some(ubv) -> min ubv.[i] v
            match x.lb with None -> c1 | Some(lbv) -> max lbv.[i] c1


        v |> Math.Vector.mapi (fun i v -> clip i v)



type stop = SmallGrad | SmallDp | MaxIters | Singular | SmallErr | InvalidObjective | NoFurther | NotDone

type lmState = {
    p : vector;
    mu : float option;
    nu : float;
    stop : stop;
    iterations : int;
}

type lmResult = {
    p : vector;
    stop : stop;
    iterations : int;
}



let epsilon = 1e-12
let isFinite v = (not (System.Double.IsInfinity v)) && (not (System.Double.IsNaN v))

let lm (f : vector -> vector) (p0 : vector)  (opts : OptimizationOptions) =

    // Starting point    
    let e0 = f p0   

    // Dimensions in objective function output
    let n = e0.Length

    // Dimensions in solution (number of function params)       
    let m = p0.Length

    // Compute initial mu
    let initMu (jacTjac : Matrix<float>) = (jacTjac.Diagonal().Maximum()) * opts.Tau

    // Take on LM step
    let stepP (state : lmState) =
        let eVect = f (state.p)
        let e = new DenseVector(eVect.InternalValues)
        let e_L2 = e.Norm(2.0)

        let mutable jac = jacobian f state.p eVect opts.Delta 
        let jacT = jac.Transpose()

        let jacTjac = jacT * jac
        let jacTe = jacT * e

        let j = ref 0
        let numActive = ref 0

        let jacTe_inf = jacTe.AbsoluteMaximum()  //Vector.fold (fun acc v -> max acc (abs v)) 0.0 jacTe
        let p_L2 = L2 state.p

        
        let mu = match state.mu with Some(mu) -> mu | None -> initMu jacTjac
        let nu = state.nu

        // Modify jacTjac
        let dia = jacTjac.Diagonal()
        for i in 0..(m-1) do
            jacTjac.[i,i] <- dia.[i] + mu
        
        let luSolve = new DenseLU(jacTjac :?> DenseMatrix)
        let Dp = -luSolve.Solve(jacTe)
        let Dp_L2 = Dp.Norm(2.0)
        let DpVect = Vector.ofArray(Dp.ToArray())

        // NOTE -- we don't know if this was solved successfully or not!  Figure out if issolved == true
        // Assume it worked for now        

        // Constrain parameter
        let p = state.p
        let pDp = p + DpVect
        

        // Compute the objective at the new location
        let eDpVect = f pDp
        let eDp = new DenseVector(eDpVect.InternalValues)
        let eDp_L2 = eDp.Norm(2.0)

        // Score improvement
        let dF = e_L2 - eDp_L2
        let dL = Dp.DotProduct(mu * Dp - jacTe)

        //printfn "e: %A" e
        //printfn "dF: %A" dF
        //printfn "dL: %A" dL 

        // Check stopping conditions
        let stopCheck =
            if e_L2 < opts.Eps3 then SmallErr
            else if jacTe_inf < opts.Eps1 then SmallGrad
            else if Dp_L2 <= opts.Eps2 * p_L2 then 
                //printfn "SmallDp. err: %A, Dp: %A" eVect Dp
                //printfn "jac: %A jac" jac
                SmallDp
            else if Dp_L2 >= (p_L2+opts.Eps2)/(epsilon**2.0) then Singular
            else if not (isFinite eDp_L2) then InvalidObjective
            else if nu > 1e8 then NoFurther
            else NotDone

        // return updates
        if dF > 0. && dL > 0. then
            { mu = Some(mu * (max 0.33 (1. - (2. * dF / dL - 1.)**3.)));
              nu = 2.0;
              p = pDp;
              stop = stopCheck;
              iterations = state.iterations + 1;
            }
          else
            { mu = Some(mu * state.nu);
              nu = state.nu * 2.0;
              p = state.p;
              stop =  stopCheck
              iterations = state.iterations + 1;
            }

    let rec lmLoop (state : lmState) k =
        if k = opts.MaxIters then 
            { state with stop = MaxIters } 
        else
            // Take a LM step and see what happened
            let r = stepP state
            match r.stop with 
                | NotDone -> lmLoop r (k+1)
                | _ -> r


    let start = { mu = None; nu = 2.0; p = p0; stop = NotDone; iterations = 0 }
    
    lmLoop start 0



let lmLinearConstraints (A : matrix) (b : vector) (f : vector -> vector) (p0 : vector)  (opts : OptimizationOptions) =
    let (Y, Z, c) = elim A b

    let bP0 = A * p0
    let d = b - bP0
    if d |> Seq.exists (fun v -> v > 1e-3)  then failwith("initial argument doesn't satisfy constraint")

    let pp = Z.Transpose * (p0 - c)

    // Find a solution the reduced objective function
    let fNew (x_z : vector) = 
        let x = c + Z * x_z
        f x

    let newResult = lm fNew pp opts

    // Convert the result back to the full solution
    { newResult with p = c + Z * newResult.p }

type ct = Interval of float * float | Low of float | High of float

let lmBoxConstrains (lb : vector option) (ub : vector option) (A : matrix) (b : vector) (f : vector -> vector) (p0 : vector)  (opts : OptimizationOptions) =
    let (Y, Z, c) = elim A b
    let n = p0.Length // Size of original problem
    let inf = System.Double.PositiveInfinity
    let w = 5e1
    
    let bP0 = A * p0
    let d = b - bP0
    if d |> Seq.exists (fun v -> v > 1e-3)  then failwith("initial argument doesn't satisfy constraint")
    
    let pp = Z.Transpose * (p0 - c)

    let constraints = Seq.init n id |> Seq.choose (fun i ->
        match lb, ub with 
            | Some(lb), Some(ub) when lb.[i] > -inf && ub.[i] < inf -> Some(Interval(lb.[i], ub.[i]), i)
            | Some(lb), _ when lb.[i] > -inf -> Some(Low(lb.[i]), i)
            | _, Some(ub) when ub.[i] < inf -> Some(High(ub.[i]), i)
            | _ -> None) |> Array.ofSeq

    let augment v boxConst =
        match boxConst with
            | Interval(low, high) -> 
                let tmp = (2.0 * v - (low + high)) / (high-low)
                w * Math.Max(tmp**2.0-1.0, 0.0)
            | Low(low) -> w * Math.Max(low - v, 0.0)
            | High(high) -> w * Math.Max(v - high, 0.0)
            
    // Find a solution the reduced objective function
    let fNew (x_z : vector) = 
        let x = c + Z * x_z
        let mainResults = (f x) |> Math.Vector.toArray

        let constraintResults = constraints |> Array.map (fun (con, idx) -> augment (x.[idx]) con)
        Array.append mainResults constraintResults |> Vector.ofArray

    let newResult = lm fNew pp opts

    // Convert the result back to the full solution
    { newResult with p = c + Z * newResult.p }