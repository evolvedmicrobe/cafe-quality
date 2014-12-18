import os
import glob
import sys


class DataofNone(object):
    def __init__(self):
            self.dic = {}
            self.blank=",".join(["NA"]*8)
    def __getitem__(self, key):
        if self.dic.has_key(key):
            return self.dic[key]
        else:            
            return self.blank
    def __setitem__(self, key, value):
        self.dic[key] = value

class FillDict(object):
    def __init__(self):
        self.dic = {}
    def __getitem__(self, key):
        if self.dic.has_key(key):
            return self.dic[key]
        else:
            new = DataofNone()
            self.dic[key]=new
            return new
        
direc = direc = sys.argv[1] # r"/Users/nigel/git/cafe-quality/master_full/"
ccs = glob.glob(direc+"*ccs.csv")
reads = glob.glob(direc+"*read_report.csv")[0]

#write header
prefix = "_".join((reads.split('/')[-1]).split("_")[:2])
header = open(reads).readline().strip() + "," + open(ccs[0]).readline().strip() + "\n"
ofn = direc+prefix+"_combined_reads.csv"
of = open(ofn,'w')
of.write(header)
          
#load error data
movies = FillDict()
d = open(reads).readlines()
for j in d[1:]:
    sp = j.split(",")
    mvn = sp[0]
    zmw = sp[1]
    movies[mvn][zmw] = j.strip()

#output all data
for read in ccs:
    data = open(read).readlines()
    for line in data[1:]:
        sp = line.split(",")
        mvn = sp[0]
        zmw = sp[1]
        data = movies[mvn][zmw]
        cmplt = data + "," + line
        of.write(cmplt)
of.close()




