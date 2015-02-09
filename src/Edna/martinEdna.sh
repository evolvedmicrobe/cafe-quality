#!/bin/bash
#LD_LIBRARY_PATH=$SEYMOUR_HOME/analysis/lib/mono/Edna:$LD_LIBRARY_PATH exec $SEYMOUR_HOME/redist/mono_fsharp/targetos/bin/mono --debug $MONO_OPTIONS $SEYMOUR_HOME/analysis/lib/mono/Edna/Edna.exe "$@"
#
# Edna from secondary.  As of May 2012, Edna is built from primary, then deployed to
# on the shared cluster FS.  A couple special symlinks control which version of mono
# and edna will be invoked by this script.  This script just sets up the environment
# and passes the arguments through
#

# clear LD_LIBRARY_PATH to ensure we use correct libs
unset LD_LIBRARY_PATH

# Scripts we use to set up the environment and run
export SCRIPT_PATH=/mnt/data3/scratch/Data/pipeline/scripts

# This is a symlink to the mono build we want to use for edna on the cluster
export MONO_PREFIX=/mnt/primary/mono/edna

. $SCRIPT_PATH/mono-profile

# This is a symlink to the edna that we want to use with smrtpipe
export EDNA_PATH=/mnt/data3/scratch/Data/edna/releases/martin
export LD_LIBRARY_PATH=$EDNA_PATH:$LD_LIBRARY_PATH

# Run edna with llvm for speed (seems to be ~20% faster?)
mono --llvm  $EDNA_PATH/Edna.exe $@

