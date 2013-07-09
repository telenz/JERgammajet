#!/bin/bash
start=$(date +%s) ## get current time/date to $start


# ----------------------------------------------------------------------------------------------------
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )
# ----------------------------------------------------------------------------------------------------
# Evaluate all systematic Uncertainties
# ----------------------------------------------------------------------------------------------------
cd scripts/
root -b -q evalJECUncertainty.C+"()"
root -b -q evalFlavorUncertainty.C+"(\"phys\")"
root -b -q evalPUUncertainty.C+"()"
root -b -q evalQCDUncertainty.C+"(1)"
rm *.C~
rm *.so
rm *.d
rm *.sh~
cd ../
# ----------------------------------------------------------------------------------------------------


# Postprocess the output
# ----------------------------------------------------------------------------------------------------
root -b -q postprocessingSysError.C+"()"
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