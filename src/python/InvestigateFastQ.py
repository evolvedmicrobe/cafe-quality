# Script to examine HP errors in the HP template
from Bio import SeqIO
import glob
fs = glob.glob("/Users/nigel/git/cafe-quality/src/PacBio.ConsensusTools/bin/Release/Test/m141008*.ccs.fastq")



def FindAcc(seqFull, bp, start, length):
    qv = seqFull.letter_annotations['phred_quality']
    seq = str(seqFull.seq)
    zmw = seqFull.id.split("/")[1]
    for j in range(start - 3, start + length + 10):
        if seq[j] == bp:
            startl = j
            end = startl + 1
            while seq[end] == bp:
                end += 1
            acc = min(qv[startl:end])
            calledLength =  end - startl
            correct = calledLength == length
            if abs(calledLength - length) < 2:
                return [bp, length, acc, correct, zmw]
            return None
    return None


forwards = [["G", 29, 5],
            ["G", 44, 10],
            ["T", 75, 10],
            ["T", 95, 5]]

reverses  = [["A", 36, 5],
             ["A", 51, 10],
             ["C", 82, 10],
             ["C", 102, 5]]


results = []
for fname in fs:
    data = list(SeqIO.parse(fname, "fastq"))
    hp = [x for x in data if abs(len(x) - 137) < 5]
    for x in hp:
        seq = str(x.seq)
        toUse = reverses
        if seq.count("TTTTTTTT") > 0:
            toUse = forwards
        for f in toUse:        
            res = FindAcc(x, f[0], f[1], f[2])
            if res:
                results.append(res)

of = open("/Users/nigel/git/cafe-quality/NotTracked/qv_compare.csv", 'w')
of.write("BP,Size,QV,Correct,ZMW\n")
for r in results:
    of.write(",".join([str(x) for x in r]) + "\n")
of.close()