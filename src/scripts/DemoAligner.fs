//#I @"C:\git\cafe-quality\lib\\"

//#r "/Users/nigel/git/cafe-quality/lib/" 
//"Bio.dll"
//"VariantCaller.dll"

open System
open Bio
open Bio.Algorithms.Alignment;
open VariantCaller;
open System.Collections.Generic
open System.Linq

[<EntryPoint>]
let main args =
    let ref = new Sequence(DnaAlphabet.Instance, "AAAAAAAAAAACCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCAAAAAAAAAAA");
    let query = new Sequence(DnaAlphabet.Instance, "AAAAAAAAAAACCCCCCCCCCCCAAAAAAAAAAA");

    let refg = new Reference(ref);
    let alns = refg.AlignSequence(query);
    let aln1 = alns.[0]
    Console.WriteLine aln1.FirstOffset
    0



