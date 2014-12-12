#!/bin/bash
start=$(date +%s) ## get current time/date to $start

ini ROOT
module load root

# ----------------------------------------------------------------------------------------------------
# Evaluate all systematic Uncertainties
# ----------------------------------------------------------------------------------------------------

cd ../
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
source runMacro.sh
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc

if [ -d plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_mc_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data_SAVED
fi
cp -r plots_2012/PF_L1CHS/data/root_files plots_2012/PF_L1CHS/systematicUncertainties/root_files_FINAL_data

cd scriptsForSystematicUncertainties/

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flavor Uncertainty Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
source runFlavorMacro.sh
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PU Uncertainty Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
source runPUMacro.sh
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JEC Uncertainty Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
source runJECMacro.sh
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QCD Uncertainty Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
source runQCDMacro.sh
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MC Uncertainty Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
source runMCMacro.sh
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
rm *.sh~
# ----------------------------------------------------------------------------------------------------


stopped=$(date +%s)
thetime=$(($stopped-$start)) ## doing some math in shell calculating time difference
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )