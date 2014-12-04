""" This is a super hacky build script designed to take a git branch and benchmark it on a set of data, 
    once it works and the problem is better defined, I'll clean it up.  Until then, here be gross code! """


import os
import subprocess
import sys

plat = sys.platform
# Hard coded paths
if plat=="darwin":
    print "Running on Mac"
    test_top_dir = "/Users/nigel/CCS_P6_C4/TestRun/"
    fofn = "/Users/nigel/CCS_P6_C4/input.fofn"
elif plat=="linux2":
    print "Running on Unix"
    fofn = "/home/UNIXHOME/ndelaney/ccswork/CCS_P6_C4/input.fofn"
    test_top_dir = "/home/UNIXHOME/ndelaney/ccswork/CCS_P6_C4/TestRun/"


# should be /Users/nigel/git/cafe-quality/src/python
start_dir = os.getcwd()
os.chdir("../")
src_dir = os.getcwd()


def GetGitInfo():
    cmd = ["git","log"]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[:3]
    op = output[0].split("\n")
    hash = op[0].split(" ")[1][:7]
    print hash
    message = op[4].strip()
    print "Building for " + hash
    print message
    return (hash, message)

def CreateTestDirectory():
    info = GetGitInfo()
    dirName = info[0]
    dirName = os.path.join(test_top_dir,dirName)
    if os.path.exists(dirName):
        os.system("rm -r "+dirName)
    os.mkdir(dirName)
    chemdir = os.path.join(dirName,"Chemistry")
    os.mkdir(chemdir)
    return dirName   

def BuildManaged(outputDir):
    cmd = "xbuild /p:Configuration=Release CafeQuality.sln"
    os.chdir(src_dir)
    res = os.system(cmd)
    if res !=0:
        raise "Failed to build managed code"
    releaseDirs = [x[0] for x in os.walk(src_dir) if x[0].count("/bin/Release")==1]
    print releaseDirs
    for j in releaseDirs:
        os.system("cp -r " + j + "/* " + test_dir + "/")
        

def MoveChemistry(test_dir):
    os.chdir(src_dir)
    mp_name = "mapping.xml"
    np = os.path.join(test_dir,"Chemistry",mp_name)
    op = "../lib/Chemistry/" + mp_name
    cmd = "cp " + op + "  " + np
    res = os.system(cmd)
    if res != 0:
        raise "Failed to move chemistry file"
    
    # Now move the training file
    ch_file = "CCSParameters.ini"
    op = "../data/" + ch_file
    np = os.path.join(test_dir,"Chemistry", ch_file)
    cmd = "cp " + op + "  " + np
    res = os.system(cmd)
    if res != 0:
        raise "Failed to move parameter file"

   
def RunTest(dir_to_run, fofn):
    cmd_base = "mono PacBio.ConsensusTools.exe"
    cmd_base += " " + fofn
    os.chdir(dir_to_run)
    res = os.system(cmd_base)
    if res != 0:
        raise "CCS Failed"


test_dir = CreateTestDirectory() 
BuildManaged(test_dir)
MoveChemistry(test_dir)
RunTest(test_dir, fofn)

