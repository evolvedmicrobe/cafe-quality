import os
direc ="/Users/nigel/git/cafe-quality/data/"
f = open("/Users/nigel/git/cafe-quality/data/corrected_ratio_ccs.fa").readlines()
o = open(direc+"corrected_ratio_ccs.trimmed.fa",'w')
c = 0
for l in f:
    if l.startswith(">"):
        c +=1
    if c < 10000:
        o.write(l)
    else:
        break
o.close()
