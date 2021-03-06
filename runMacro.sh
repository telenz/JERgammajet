#!/bin/bash
start=$(date +%s) ## get current time/date to $start


# ----------------------------------------------------------------------------------------------------
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
# ----------------------------------------------------------------------------------------------------

# Run whole analysis for data
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_data.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
# ----------------------------------------------------------------------------------------------------

# Postprocess the output
# ----------------------------------------------------------------------------------------------------
root -b -l scripts << EOF
.L DataMCComparison.C+
plotDataMCComparisonFINAL()
EOF
# ----------------------------------------------------------------------------------------------------
# Remove all unnecessary files
rm *.C~
rm CODE/*.C~
rm CODE/*.h~
rm scripts/*.C~
rm scripts/*.h~
rm *.so
rm *.d
rm scripts/*.so
rm scripts/*.d
rm *.sh~
# ----------------------------------------------------------------------------------------------------


stopped=$(date +%s)
thetime=$(($stopped-$start)) ## doing some math in shell calculating time difference
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )