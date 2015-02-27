module Caller16S
(* Script to take 16S sequences and generate alignments to show where we need to mask them. *)
open System
open System.IO
open VariantCaller
open Bio
open Bio.IO.FastA
open Bio.Algorithms.Alignment
open PacBio.Data
open PacBio.Utils
open ConsensusCore
open LoadZMWs
open System.Linq
open PacBio.Consensus
open System.Collections.Generic
open Microsoft.FSharp.Collections


let seqs = 
   let fr = new Bio.IO.FastA.FastAParser("/Users/nigel/git/cafe-quality/NotTracked/unique_16s.fasta")
   fr.Parse().ToList()

let getName (id:ISequence) = id.ID.Split('.').[1]

let processGroup (grp : string * seq<ISequence>) =
    Console.WriteLine (fst grp)
    let sp = fst grp
    let s = (snd grp).ToList();
    let ref = new Reference(s.[0] :?> Sequence)
    let alnAndCall (s2 : ISequence) = 
        let aln = ref.AlignSequence(s2 :?> Sequence)
        Console.WriteLine(aln.[0])
        let variants = VariantCaller.VariantCaller.CallVariants(aln.[0], ref.RefSeq);
        Console.WriteLine("Variants = " + variants.Count.ToString())
        variants |> Seq.map (fun v -> Console.WriteLine(v.ToString())) |> Seq.length
        1
    s |> Seq.skip 1 |> Seq.map alnAndCall |> Seq.sum


let main =
   let groups = Seq.groupBy getName seqs
   groups |> Seq.map processGroup |> Seq.length
   0
