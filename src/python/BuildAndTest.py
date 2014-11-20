""" This is a super hacky build script designed to take a git branch and benchmark it on a set of data, 
    once it works and the problem is better defined, I'll clean it up.  Until then, here be gross code! """


import os
import subprocess

# Hard coded paths
test_top_dir = "/Users/nigel/CCS_P6_C4/TestRun/"
fofn = "/Users/nigel/CCS_P6_C4/input.fofn"

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
        
def BuildUnmanaged():
    cc_dir = os.path.join(src_dir, "ConsensusCore")
    os.chdir(cc_dir)
    res = os.system("make")
    if res != 0:
        raise "Failed to build unmanaged code"

def RunTest(dir_to_run, fofn):
    cmd_base = "mono PacBio.ConsensusTools.exe CircularConsensus -n 8"
    outName = dir_to_run.split("/")[-1]
    cmd_base += " -o " + outName
    cmd_base += " -fofn=" + fofn
    os.chdir(dir_to_run)
    res = os.system(cmd_base)
    if res != 0:
        raise "CCS Failed"


test_dir = CreateTestDirectory() 
BuildUnmanaged()
BuildManaged(test_dir)
RunTest(test_dir, fofn)

