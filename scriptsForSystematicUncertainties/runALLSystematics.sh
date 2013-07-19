#!/bin/bash
start=$(date +%s) ## get current time/date to $start

# ----------------------------------------------------------------------------------------------------
# Evaluate all systematic Uncertainties
# ----------------------------------------------------------------------------------------------------
source runFlavorMacro.sh
source runPUMacro.sh
source runJECMacro.sh
source runQCDMacro.sh
source runMCMacro.sh
rm *.sh~
# ----------------------------------------------------------------------------------------------------


stopped=$(date +%s)
thetime=$(($stopped-$start)) ## doing some math in shell calculating time difference
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )