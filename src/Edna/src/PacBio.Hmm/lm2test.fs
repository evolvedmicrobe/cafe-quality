module lm2test

open System
open System.Diagnostics
open System.Reflection
open System.IO
open System.Text    

open NUnit.Framework
open NUnit.Framework.Constraints
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics
open lm2

// Test objective functions for various flavors of LM optimization.

// Rosenbrock function, global minimum at (1, 1)
let ros (p : Vector<float>) =
    let ROSD = 25.0
    let v = ((1.0-p.[0])*(1.0-p.[0]) + ROSD*(p.[1]-p.[0]*p.[0])*(p.[1]-p.[0]*p.[0]));
    vector [5.0; v]

// Modified Rosenbrock - global minimum at (1,1)
let modros (p : Vector<float>) =
    let modroslam = 1e2
    let x = vector [3.0; 0.0]
    x.[0] <- 10.0 * (p.[1] - p.[0] * p.[0])
    x.[1] <- 1.0 - p.[0]
    x.[2] <- modroslam

    x

let powell (p : Vector<float>) =
    let x = vector [3.0; 0.0]
    x.[0] <- p.[0]
    x.[1] <- 10.0 * p.[0] / (p.[0] + 0.1) + 2. * p.[1] * p.[1]
    x

// Wood's function, minimum at (1, 1, 1, 1)
let wood (p : Vector<float>) =
    let x = vector [6.0; 0.0]
    x.[0] <- 10. * (p.[1] - p.[0] * p.[0])
    x.[1] <- 1.0 - p.[0]
    x.[2] <- sqrt(90.) * (p.[3] - p.[2] * p.[2])
    x.[3] <- 1.0 - p.[2]
    x.[4] <- sqrt(10.0) * (p.[1] + p.[3] - 2.0)
    x.[5] <- (p.[1] - p.[3]) / sqrt(10.0)
    x


// Meyer's (reformulated) problem, minimum at (2.48, 6.18, 3.45) */
let meyer (p : Vector<float>) = 
    let n = 16;
  
    let xx = [| 
        34.780; 28.610; 23.650; 19.630; 16.370; 13.720; 11.540; 9.744;
        8.261; 7.030; 6.005; 5.147; 4.427; 3.820; 3.307; 2.872 |]
    DenseVector.init xx.Length (fun i ->
        let ui = 0.45 * 0.05 * (float i)
        p.[0] * Math.Exp(10.0 * p.[1] / (ui + p.[2]) - 13.0) - (float xx.[i]))


[<TestFixture>]
type TestUtils() =

    [<Test>]
    member x.Rosenbrock() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = vector [2.0; 0.5] 
        let result = lm2.lm ros start opts
        printfn "%A" result


    [<Test>]
    member x.ModifiedRosenbrock() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = vector [2.0; 0.5] 
        let result = lm2.lm modros start opts

        printfn "%A" result

    [<Test>]
    member x.Powell() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = vector [2.0; 0.5] 
        let result = lm2.lm powell start opts

        printfn "%A" result

    [<Test>]
    member x.Wood() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = vector [4.0; 0.1] 
        let result = lm2.lm wood start opts
        printfn "%A" result


    [<Test>]
    member x.Meyer() =
        let opts = lm2.OptimizationOptions.standard(1000)
        let start = vector [3.0; 0.1] 
        let result = lm2.lm meyer start opts

        printfn "%A" result
  