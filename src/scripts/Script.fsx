File to count the length of sequences seen by CCS

#I @"C:\git\cafe-quality\lib\"
#r @"Bio.dll"
#r @"RDotNet.dll"
#r @"RDotNet.FSharp.dll"
#r @"RDotNet.NativeLibrary.dll"
#r @"RProvider.dll"


open System.IO
open System
open System.Linq
open RDotNet
open RProvider
open RProvider.``base``
open RProvider.graphics
open RProvider.stats

open Bio
open Bio.IO.FastA

let dataDir = @"C:\git\cafe-quality\data"

let files = DirectoryInfo(dataDir).GetFiles()

let zf = files |> Seq.where (fun z -> z.Name.EndsWith("ccs.fasta.gz")) |> Seq.toList

let getSeqs (fname: FileInfo) = 
   let reader = new FastAZippedParser(fname.FullName)
   reader.Parse() |> Seq.toList

let points = zf |> Seq.collect getSeqs |> Seq.toArray
Define your library scripting code here

let sizes = points |> Seq.map (fun h -> Math.Log10((float)h.Count)) |> Seq.toArray

let sizes = points |> Seq.map (fun h -> (float)h.Count) |> Seq.toArray

let mkPlot data = 
   let args = new System.Collections.Generic.Dictionary<string,Object>()
   args.["x"] <- data
   args.["main"] <-"Log 10 size distribution"
   args.["col"] <-"blue"
   args.["breaks"]<-200
   args.["xlab"]<-"Size "
   R.hist(args)
mkPlot sizes

let largest = points |> Seq.maxBy (fun u -> u.Count)
(largest :?> Sequence).ConvertToString()