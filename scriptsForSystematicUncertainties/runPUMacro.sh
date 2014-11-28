#!/bin/bash
start=$(date +%s) ## get current time/date to $start

cd ../

# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the PU uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) Deactivate trigger requirement
sed -i 's/applyTriggeronMC=true;/applyTriggeronMC=false;/g' CODE/myDeclarations.h
sed -i 's/PUreweighting=false;/PUreweighting=true;/g' CODE/myDeclarations.h
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
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue

# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (upward variation)
# ----------------------------------------------------------------------------------------------------
sed -i 's/MBXS69400/MBXS73000/g' CODE/readGammaJet.C
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (downward variation)
# ----------------------------------------------------------------------------------------------------
sed -i 's/MBXS73000/MBXS65800/g' CODE/readGammaJet.C
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation
# ----------------------------------------------------------------------------------------------------

# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) Deactivate trigger requirement
sed -i 's/applyTriggeronMC=false;/applyTriggeronMC=true;/g' CODE/myDeclarations.h
# 2) Correct Pu weights
sed -i  's/MBXS65800/MBXS69400/g' CODE/readGammaJet.C

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

cd scriptsForSystematicUncertainties/

stopped=$(date +%s)
thetime=$(($stopped-$start)) ## doing some math in shell calculating time difference
echo "Time needed: " $(date -d "1970-01-01 $thetime sec" +"%H:%M:%S") / $thetime "secs" ## prints out the time ( like: 00:01:29 / 89 secs )