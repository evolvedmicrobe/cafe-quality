module lm2
(*  Module to compute levenberg-marquadt optimization of a given function.
    Note from Nigel: Based on variable names and the unit tests, this 
    appears to be an F# port of the levmar library, in particular lm_core.c

    *References*
    Wiki: https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
    LevMar: http://users.ics.forth.gr/~lourakis/levmar/
    LevMar Algorithm Description: http://users.ics.forth.gr/~lourakis/levmar/levmar.pdf
    Recent Review of LM Options: http://arxiv.org/abs/1201.5885v1
*)

open MathNet.Numerics.LinearAlgebra.Double.Solvers
open MathNet.Numerics.LinearAlgebra.Double
open MathNet.Numerics.LinearAlgebra
open System
open PacBio.Hmm.Utils

let L2 = Vector.fold (fun acc v -> acc + v**2.) 0.0
type OptimizationOptions = {
    (* Lowest step size used for computing numeric derivatives 
     when calculating the jacobian. *)
    Delta : float;
    (* EPS Values are stopping conditions. See description in levmar. *)        
    Eps1 : float; //  ||J^T e||_inf the infinity norm, e.g. highest absolute value 
    Eps2 : float; // ||Dp||_2 is the amount by which p is currently being shifted at each iteration; 
    Eps3 : float; //||e||_2 is a measure of the error between the model function at the current estimate for p and the data.
    // Scale factor for initial dampening, set by levmar
    Tau : float;
    // Maximum LM iterations allowed
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

MathNet.Numerics.Control.UseSingleThread()


(* Compute a numerically calculated jacobian matrix using a delta of either the parameter given, or 1e-4 the current parameter value *)
//TODO: Avoid numeric derivatives?
let jacobian (f : Vector<float> -> Vector<float>) (p : Vector<float>) (hx : Vector<float>) delta  =
    let jac = DenseMatrix.zero<float> hx.Count p.Count
    for j in 0 .. (p.Count - 1) do
        let d = max delta (p.[j] * 1e-4)
        let tmp = p.[j]
        p.[j] <- p.[j] + d
        let hxx = f p
        p.[j] <- tmp
        for i in 0 .. (hx.Count - 1) do
            jac.[i,j] <- (hxx.[i] - hx.[i]) / d // dydx
    jac


type stop = SmallGrad | SmallDp | MaxIters | Singular | SmallErr | InvalidObjective | NoFurther | NotDone

type lmState = {
    p : Vector<float>;
    mu : float option;
    nu : float;
    stop : stop;
    iterations : int;
}

type lmResult = {
    p : Vector;
    stop : stop;
    iterations : int;
}

let epsilon = 1e-12
let isFinite v = (not (System.Double.IsInfinity v)) && (not (System.Double.IsNaN v))

(* Perform a Levenberg-Marquadt minimization of the following function using
the given error function, parameter vector and optimization options.
Note: The error should directly calculate the residuals (unsquared), not the sum of squares, the 
squared errors or the predicted values.

Note: If profiling shows this as an Edna bottleneck, lots of possible improvements, including:
    1 - Avoid recomputing matrices after rejected steps
    2 - Modify arrays in place rather than copy.
    3 - Use Cholesky instead of LU (possible accuracy loss, but we are using floats).
 *)
let lm (f : Vector<float> -> Vector<float>) (p0 : Vector<float>)  (opts : OptimizationOptions) =

    (* Compute initial dampening factor,
       setting it near the maximum of the diagonal gives a rough balance 
       between steepest descent and gauss-newton *)
    let initMu (jacTjac : Matrix<float>) = (jacTjac.Diagonal().Maximum()) * opts.Tau

    // Take one LM step
    let stepP (state : lmState) =
        let error = f (state.p) //compute current error 
        let error_L2 = error.L2Norm() //compute penalty

        // Numerically compute jacobian
        let jac = jacobian f state.p error opts.Delta 

        let jacT = jac.Transpose()
        //TODO: Consider speeding up multiplication as this is symmetric? 
        let jacTjac = jacT * jac  //Hessian approximation
        let jacTe = jacT * error  

        let j = ref 0
        let numActive = ref 0
        //value to check as stopping condition 
        let jacTe_inf = jacTe.AbsoluteMaximum()  
        let p_L2 = L2 state.p // Get sum of squared parameters, equivalent to the norm without the square root
                
        // Initialize the dampening factor
        let mu = match state.mu with Some(mu) -> mu | None -> initMu jacTjac
        let nu = state.nu

        // Modify jacTjac to augment normal equations
        let dia = jacTjac.Diagonal()
        for i in 0..(jacTjac.RowCount-1) do
            (*  FIXME: 
                As written this appears to not use the Marquadt dampening matrix.
                This decision propogates from the original levmar library and in turn the lecture notes by 
                K. Madsen, H.B. Nielsen, O. Tingleff levmar is based on.  In particular, it seems the program
                attempts to calculate the update by solving:
                (jacTjac + \mu I)Dp = jacTe rather than:
                (jacTjac + \mu diag(jacTjac)DP = jacTe as is described on wikipedia.
                The later framework (which involves a multiplication rather than a simple addition) is used by 
                numerical recipes and some other routines.  However, this is still a valid algorithm and appears
                related to the approach in the old fortran routine minpack, with some justification for it (and a
                possible improvement) given in this 2012 article: http://arxiv.org/abs/1201.5885v1 
                (see section 2.2 for a discussion of choosing the dampening matrix). 
                Not altering at present as it appears to work...
            *)
            jacTjac.[i,i] <- dia.[i] + mu
        
        // Solve for parameter update, (jacTjac + \mu I)Dp = jacTe, using LU decompostion
        //FIXME: We don't know if this was solved successfully or not!  Figure out if issolved == true
        let Dp = - (jacTjac.LU().Solve(jacTe))
        let Dp_L2 = Dp.L2Norm()

        // update parameter
        let pDp = state.p + Dp        

        // Compute the objective at the new location
        let new_error = f pDp
        let new_error_L2 = new_error.L2Norm()

        // Score improvement
        let dF = error_L2 - new_error_L2
        let dL = Dp.DotProduct(mu * Dp - jacTe)

        //printfn "error: %A" error
        //printfn "dF: %A" dF
        //printfn "dL: %A" dL 

        // Check stopping conditions
        let stopCheck =
            if error_L2 < opts.Eps3 then SmallErr
            else if jacTe_inf < opts.Eps1 then SmallGrad
            else if Dp_L2 <= opts.Eps2 * p_L2 then 
                //printfn "SmallDp. err: %A, Dp: %A" eVect Dp
                //printfn "jac: %A jac" jac
                SmallDp
            else if Dp_L2 >= (p_L2+opts.Eps2)/(epsilon**2.0) then Singular
            else if not (isFinite new_error_L2) then InvalidObjective
            else if nu > 1e8 then NoFurther
            else NotDone

        // reduction in error, increment is accepted
        if dF > 0. && dL > 0. then
            { mu = Some(mu * (max 0.33 (1. - (2. * dF / dL - 1.)**3.)));
              nu = 2.0;
              p = pDp;
              stop = stopCheck;
              iterations = state.iterations + 1;
            }
          else //reject step, note we could avoid recomputing jacobian, hessian, etc. after this.
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


    let start = { mu = None; nu = 2.0; p = p0; stop = NotDone; iterations = 0;}
    
    lmLoop start 0