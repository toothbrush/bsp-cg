#!/bin/bash

CMD="$@"
LOGNAME=`date +'%s'`
LOGNAME="$LOGNAME-$$"

# redirect all output following this line to logfiles.

exec 1> $LOGNAME.out.log 2>$LOGNAME.err.log
echo $CMD
sh -c "time $CMD"
