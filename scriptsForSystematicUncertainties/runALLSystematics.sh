#!/bin/bash
start=$(date +%s) ## get current time/date to $start

# ----------------------------------------------------------------------------------------------------
# Evaluate all systematic Uncertainties
# ----------------------------------------------------------------------------------------------------

cd ../
source runMacro.sh
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc

if [ -d plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data_SAVED
fi
cp -r plots_2012/PF_L1CHS/data/root_files plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data

cd scriptsForSystematicUncertainties/

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