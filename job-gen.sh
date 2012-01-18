#!/bin/bash

if [ "$1" = "" ] ; then
echo Provide an argument.
fi

mats=$1
for i in $mats/*.emm ; do

for p in 1 2 4 8 ; do

mat="$i-P$p"
u="$i-u$p"
v="$i-v$p"

jobname=`basename $mat`.job

tee $jobname > /dev/null <<HERE
# @ node = 1
#
# @ tasks_per_node = $p
# @ notification = never
# @ input = /dev/null
# @ output = out-nz-$mat.\$(jobid)
# @ error = err-nz-$mat.\$(jobid)
# @ wall_clock_limit = 00:30:00
# @ job_type = parallel
#
# @ queue
#
# @ node_usage = shared
cd \$HOME/Students10/pdwalt/bsp-cg/
./src/cg $mat $u $v
HERE


done
done
