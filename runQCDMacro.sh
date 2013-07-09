#!/bin/bash
start=$(date +%s) ## get current time/date to $start


# ----------------------------------------------------------------------------------------------------
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/addQCD=false;/addQCD=true;/g' "$i"
root -b -q jetphoton_mc.C+"(10**9,1)"
root -b -q jetphoton_mc.C+"(10**9,2)"
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/mc/root_files_QCDUncertainty/root_files_test5
for i in CODE/myDeclarations.h ; sed -i 's/addQCD=true;/addQCD=false;/g' "$i"
root -b -q jetphoton_mc.C+"(10**9,1)"
root -b -q jetphoton_mc.C+"(10**9,2)"
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/mc/root_files_WithoutTriggerWithPUWeightEq1_EtaBin1p3_softBinning_Neq0_otherAlphaBinning
# ----------------------------------------------------------------------------------------------------

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