module lm2test

open System
open System.Diagnostics
open System.Reflection
open System.IO
open System.Text    

open NUnit.Framework
open NUnit.Framework.Constraints

open Microsoft.FSharp.Math


open lm2

// Test objective functions for various flavors of LM optimization.

// Rosenbrock function, global minimum at (1, 1)
let ros (p : vector) =
    let ROSD = 25.0
    let v = ((1.0-p.[0])*(1.0-p.[0]) + ROSD*(p.[1]-p.[0]*p.[0])*(p.[1]-p.[0]*p.[0]));
    Vector.create 5 v

// Modified Rosenbrock - global minimum at (1,1)
let modros (p : vector) =
    let modroslam = 1e2
    
    let x = Vector.create 3 0.0
    x.[0] <- 10.0 * (p.[1] - p.[0] * p.[0])
    x.[1] <- 1.0 - p.[0]
    x.[2] <- modroslam

    x

let powell (p : vector) =
    let x = Vector.create 3 0.0
    x.[0] <- p.[0]
    x.[1] <- 10.0 * p.[0] / (p.[0] + 0.1) + 2. * p.[1] * p.[1]
    x

// Wood's function, minimum at (1, 1, 1, 1)
let wood (p : vector) =
    let x = Vector.create 6 0.0
    x.[0] <- 10. * (p.[1] - p.[0] * p.[0])
    x.[1] <- 1.0 - p.[0]
    x.[2] <- sqrt(90.) * (p.[3] - p.[2] * p.[2])
    x.[3] <- 1.0 - p.[2]
    x.[4] <- sqrt(10.0) * (p.[1] + p.[3] - 2.0)
    x.[5] <- (p.[1] - p.[3]) / sqrt(10.0)
    x


// Meyer's (reformulated) problem, minimum at (2.48, 6.18, 3.45) */
let meyer (p : vector) = 
    let n = 16;

    let xx = [| 
        34.780; 28.610; 23.650; 19.630; 16.370; 13.720; 11.540; 9.744;
        8.261; 7.030; 6.005; 5.147; 4.427; 3.820; 3.307; 2.872 |]

    Vector.init xx.Length (fun i ->
        let ui = 0.45 * 0.05 * (float i)
        p.[0] * Math.Exp(10.0 * p.[1] / (ui + p.[2]) - 13.0) - (float xx.[i]))


[<TestFixture>]
type TestUtils() =

    [<Test>]
    member x.Rosenbrock() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = Vector.create 2 0.5 
        let result = lm2.lm ros start opts

        printfn "%A" result


    [<Test>]
    member x.ModifiedRosenbrock() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = Vector.create 2 0.5 
        let result = lm2.lm modros start opts

        printfn "%A" result

    [<Test>]
    member x.Powell() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = Vector.create 2 0.5 
        let result = lm2.lm powell start opts

        printfn "%A" result

    [<Test>]
    member x.Wood() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = Vector.create 4 0.1 
        let result = lm2.lm wood start opts

        printfn "%A" result


    [<Test>]
    member x.Meyer() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = Vector.create 3 0.1 
        let result = lm2.lm meyer start opts

        printfn "%A" result


    // Hock - Schittkowski problem 28
    [<Test>]
    member x.HockSchittkowski28() =

        let m=3
        let n=3

        let f (p: vector) =
            let t1 = p.[0] + p.[1]
            let t2 = p.[1] + p.[2]
            Vector.create n (t1*t1 + t2*t2)
            
        let p = Vector.create 3 0.0
        p.[0] <- -4.0 
        p.[1] <- 1.0
        p.[2] <- 1.0

        let A = [| 1.0; 2.0; 3.0 |] |> RowVector.ofArray |> Matrix.ofRowVector
        let b = 1.0 |> Vector.ofScalar

        let opts = lm2.OptimizationOptions.standard(1000)
        let pp = lm2.lmLinearConstraints A b f p opts


        printfn "Final: %A" pp

        let target = [| 0.5; -0.5; 0.5 |]  |> Vector.ofArray
        ()


    // Hock - Schittkowski (modified) problem 52 (box/linearly constrained), minimum at (-0.09, 0.03, 0.25, -0.19, 0.03)
    [<Test>]
    member x.HockSchittkowski52() =
        let c1 = [| 1.0; 3.0; 0.0; 0.0; 0.0 |] 
        let c2 = [| 0.0; 0.0; 1.0; 1.0; -2.0 |]
        let c3 = [| 0.0; 1.0; 0.0; 0.0; -1.0 |]

        let b = Vector.create 3 0.0
        let A = Matrix.ofSeq [c1; c2; c3]


        let inf = System.Double.PositiveInfinity
        let ub = [| inf; 0.3; 0.25; 0.3; 0.3 |] |> Vector.ofArray |> Some
        let lb = [| -0.09; 0.0; -inf; -0.2; 0.0 |] |> Vector.ofArray |> Some

        let p0 = Vector.create 5 2.0

        let f (p : vector) =
            let v = Vector.create 4 0.0
            v.[0] <- 4.0 * p.[0] - p.[1]
            v.[1] <- p.[1] + p.[2] - 2.0;
            v.[2] <- p.[3] - 1.0;
            v.[3] <- p.[4] - 1.0;
            v

        let opts = lm2.OptimizationOptions.standard(1000)
        let res = lm2.lmBoxConstrains lb ub A b f p0 opts

        printfn "%A" res

(* Hock - Schittkowski (modified) problem 52 (box/linearly constrained), minimum at (-0.09, 0.03, 0.25, -0.19, 0.03)
 * constr1: p[0] + 3*p[1] = 0;
 * constr2: p[2] +   p[3] - 2*p[4] = 0;
 * constr3: p[1] -   p[4] = 0;
 *
 * To the above 3 constraints, we add the following 5:
 * constr4: -0.09 <= p[0];
 * constr5:   0.0 <= p[1] <= 0.3;
 * constr6:          p[2] <= 0.25;
 * constr7:  -0.2 <= p[3] <= 0.3;
 * constr8:   0.0 <= p[4] <= 0.3;
 *
 */
void modhs52(double *p, double *x, int m, int n, void *data)
{
  x[0]=4.0*p[0]-p[1];
  x[1]=p[1]+p[2]-2.0;
  x[2]=p[3]-1.0;
  x[3]=p[4]-1.0;
}
*)
(*
/* Schittkowski (modified) problem 235 (box/linearly constrained), minimum at (-1.725, 2.9, 0.725)
 * constr1: p[0] + p[2] = -1.0;
 *
 * To the above constraint, we add the following 2:
 * constr2: p[1] - 4*p[2] = 0;
 * constr3: 0.1 <= p[1] <= 2.9;
 * constr4: 0.7 <= p[2];
 *
 */

   x[0]=0.1*(p[0]-1.0);
  x[1]=p[1]-p[0]*p[0];

 *)
    // /* Schittkowski modified problem 235 */
    [<Test>]
    member x.HockSchittkowski235() =
        let c1 = [| 1.0; 0.0; 1.0;|] 

        let b = Vector.create 1 0.0
        b.[0] <- -1.0
        let A = Matrix.ofSeq [c1]

        let inf = System.Double.PositiveInfinity
        let ub = [| inf; 2.9; inf |] |> Vector.ofArray |> Some
        let lb = [| -inf; 0.1; 0.7 |] |> Vector.ofArray |> Some

        let p0 = Vector.ofArray [| -2.0; 3.0; 1.0 |]

        let f (p : vector) =
            let v = Vector.create 2 0.0
            v.[0] <- 0.1 * (p.[0] - 1.0)
            v.[1] <- p.[1] - p.[0]*p.[0]
            v

        let opts = lm2.OptimizationOptions.standard(1000)
        let res = lm2.lmBoxConstrains lb ub A b f p0 opts

        printfn "%A" res



    [<Test>]
    member x.Qr() =
        let r = new Random(42)

        let m = 10 // number of parameters
        let n = 3 // number of constraints
        let k = 3 // number of points per constraint

        let randomConstraint () =
            let ca = Array.create m 0.0
            for i in 0 .. k do
                ca.[r.Next(m)] <- 1.0
            done
            ca

        let A = (Seq.init n (fun i -> randomConstraint())) |> Matrix.ofSeq
        let b = Vector.create n 1.0

        let (Y, Z, c) = lm2.elim A b

        // Now x = Y*b + Z*x_z = c * Z*x_z
        // and Ax=b for all x_z

        let x_z = Vector.init (m-n) (fun i -> r.NextDouble())

        let x = c + Z * x_z

        let bPrime = A*x

        printfn "b: %A" b
        printfn "b': %A" bPrime
