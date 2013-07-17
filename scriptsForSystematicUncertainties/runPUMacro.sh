#!/bin/bash
start=$(date +%s) ## get current time/date to $start

cd ../

# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the PU uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) Deactivate trigger requirement
for i in CODE/myDeclarations.h ; sed -i 's/applyTriggeronMC=true;/applyTriggeronMC=false;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/PUreweighting=false;/PUreweighting=true;/g' "$i"
# ----------------------------------------------------------------------------------------------------

# Run whole analysis for MC
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/centralValue

# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (upward variation)
# ----------------------------------------------------------------------------------------------------
for i in CODE/readGammaJet.C ; sed -i 's/MBXS69400/MBXS73000/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/upwardVariation
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (downward variation)
# ----------------------------------------------------------------------------------------------------
for i in CODE/readGammaJet.C ; sed -i 's/MBXS73000/MBXS65800/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_PUUncertainty/downwardVariation
# ----------------------------------------------------------------------------------------------------

# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) Deactivate trigger requirement
for i in CODE/myDeclarations.h ; sed -i 's/applyTriggeronMC=false;/applyTriggeronMC=true;/g' "$i"
# 2) Correct Pu weights
for i in CODE/readGammaJet.C ; sed -i  's/MBXS65800/MBXS69400/g' "$i"

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