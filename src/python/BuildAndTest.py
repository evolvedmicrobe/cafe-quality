""" This is a super hacky build script designed to take a git branch and benchmark it on a set of data, 
    once it works and the problem is better defined, I'll clean it up.  Until then, here be gross code! """


import os
import subprocess
import sys

def VerifyLinux():
    """ Silly method to remind me what to do in case of errors """
    loaded = os.environ['LOADEDMODULES']
    if loaded.count("swig") == 0:
        raise " You must load swig first"
    if loaded.count("mono") ==0:
        print "*"*20
        print "Mono not loaded by module, you know what you're at right?"
        print "*"*20
    if loaded.count("hdf5-tools/1.8.14") ==0:
        raise "You need to Load HDF5 1.8.14"


plat = sys.platform
# Hard coded paths
if plat=="darwin":
    print "Running on Mac"
    test_top_dir = "/Users/nigel/CCS_P6_C4/TestRun/"
    fofn = "/Users/nigel/CCS_P6_C4/input.fofn"
    mono_opts = ""
    src_top_dir = "/Users/nigel/git/cafe-quality/src/"
elif plat=="linux2":
    VerifyLinux()
    print "Running on Unix"
    fofn = "/home/UNIXHOME/ndelaney/ccswork/CCS_P6_C4/input.fofn"
    test_top_dir = "/home/UNIXHOME/ndelaney/ccswork/CCS_P6_C4/TestRun/"
    mono_opts = " --gc=boehm "
    src_top_dir = "/home/UNIXHOME/ndelaney/git/cafe-quality/src/"

# should be /Users/nigel/git/cafe-quality/src/python
start_dir = os.getcwd()
os.chdir("../")
src_dir = os.getcwd()

def RemoveDependencies():
    print "Removing Dependencies"
    cmd_top = "find " + src_top_dir + " | grep '\.exe'  "
    os.system(cmd_top)
    cmd = cmd_top + " | xargs rm"
    res = os.system(cmd)
    if res != 0:
        raise "Failed to remove old executables!"
    cmd_top = "find " + src_top_dir + " | grep '\.dll' | xargs rm  "
    res = os.system(cmd_top)
    if res != 0:
        raise "Could not remove old dlls!"

def GetGitInfo():
    cmd = ["git","log"]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[:3]
    op = output[0].split("\n")
    hash = op[0].split(" ")[1][:7]
    branch = GetGitBranchName()
    hash = branch + "_" + hash
    message = op[4].strip()
    print "Building for " + hash
    print message
    return (hash, message)

def GetGitBranchName():
    cmd = ["git", "branch"]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE ).communicate()[:3]
    op = [x for x in output[0].split("\n") if x.count("*") == 1]
    print op
    branch = op[0].split(" ")[1]
    return branch

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
        
def BuildUnmanaged():
    cc_dir = os.path.join(src_dir, "ConsensusCore")
    os.chdir(cc_dir)
    res = os.system("./configure")
    if res !=0:
        raise "Failed to run unmanaged configure"
    res = os.system("make")
    if res != 0:
        raise "Failed to build unmanaged code"

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

def MoveCore(test_dir):
    os.chdir(src_dir)
    mp_name = "libConsensusCore.so"
    np = os.path.join(test_dir,mp_name)
    op = "../lib/" + mp_name
    cmd = "cp " + op + "  " + np
    res = os.system(cmd)
    if res != 0:
        raise "Failed to move ConsensusCore  file"
    
def RunTest(dir_to_run, fofn):
    cmd_base = "mono " + mono_opts +  " PacBio.ConsensusTools.exe CircularConsensus -n 8 -directional -csv"
    outName = dir_to_run.split("/")[-1]
    cmd_base += " -o " + outName
    cmd_base += " -fofn=" + fofn
    os.chdir(dir_to_run)
    res = os.system(cmd_base)
    if res != 0:
        raise "CCS Failed"


test_dir = CreateTestDirectory()

RemoveDependencies()
BuildUnmanaged()
BuildManaged(test_dir)
MoveChemistry(test_dir)
MoveCore(test_dir)
RunTest(test_dir, fofn)

