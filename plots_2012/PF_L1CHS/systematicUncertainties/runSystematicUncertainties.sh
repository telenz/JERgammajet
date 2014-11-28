#!/bin/bash
start=$(date +%s) ## get current time/date to $start

ini ROOT
module load root

# ----------------------------------------------------------------------------------------------------
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )
# ----------------------------------------------------------------------------------------------------
# Evaluate all systematic Uncertainties
# ----------------------------------------------------------------------------------------------------
cd scripts/
root -b -l -q evalJECUncertainty.C+"()"
root -b -l -q evalFlavorUncertainty.C+"(\"algo\")"
root -b -l -q evalPUUncertainty.C+"()"
root -b -l -q evalQCDUncertainty.C+"(1)"
root -b -l -q evalMCUncertainty.C+"()"
rm *.C~
rm *.so
rm *.d
rm *.sh~
cd ../
# ----------------------------------------------------------------------------------------------------


# Postprocess the output
# ----------------------------------------------------------------------------------------------------
root -b -l -q postprocessingSysError.C+"()"
rm *.C~
rm *.so
rm *.d
rm *.sh~
# ----------------------------------------------------------------------------------------------------
# Remove all unnecessary files

# ----------------------------------------------------------------------------------------------------


stopped=$(date +%s)
thetime=$(($stopped-$start)) ## doing some math in shell calculating time difference
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )