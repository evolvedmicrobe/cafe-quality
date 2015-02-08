#!/bin/bash
. /mnt/software/Modules/current/init/bash
module load edna/LTS
Edna.sh -n 8 --maxTraces 1000000 $@ && exit 0

# If Edna didn't end with a success, run stuff below:
#    - first try rerunning it if this is a core dump
/bin/ls | grep core\. > /dev/null && \
    Edna.sh -n 8 --maxTraces 1000000 $@ && exit 0
#    - crash otherwise
exit $?
