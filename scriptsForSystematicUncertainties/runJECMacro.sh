#!/bin/bash
start=$(date +%s) ## get current time/date to $start

cd ../

# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the QCD uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) deactivate trigger requirement
sed -i 's/applyTriggeronMC=true;/applyTriggeronMC=false;/g' "$i"
# 2) deactivate PU weighting
sed -i 's/PUreweighting=true;/PUreweighting=false;/g' "$i"
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
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/centralValue ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/centralValue plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/centralValue_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/centralValue
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (upwards variation)
# ----------------------------------------------------------------------------------------------------
sed -i 's/corrJets.add(j,\/\*(1.+jetUncert\[j\])\*\/jetCorrL1\[j\]/corrJets.add(j,(1.+jetUncert\[j\])\*jetCorrL1\[j\]/g' CODE/applyCuts.C
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/upwardVariation ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/upwardVariation plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/upwardVariation_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/upwardVariation
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (downward variation)
# ----------------------------------------------------------------------------------------------------
sed -i 's/corrJets.add(j,(1.+jetUncert\[j\])\*jetCorrL1\[j\]/corrJets.add(j,(1.-jetUncert\[j\])\*jetCorrL1\[j\]/g' CODE/applyCuts.C
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/downwardVariation ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/downwardVariation plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/downwardVariation_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_JECUncertainty/downwardVariation
# ----------------------------------------------------------------------------------------------------

# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) deactivate trigger requirement
sed -i 's/applyTriggeronMC=false;/applyTriggeronMC=true;/g' CODE/myDeclarations.h
# 2) deactivate PU weighting
sed -i 's/PUreweighting=false;/PUreweighting=true;/g' CODE/myDeclarations.h
# 1) Switch to ak7PFCHS jets
sed -i 's/corrJets.add(j,(1.-jetUncert\[j\])\*jetCorrL1\[j\]/corrJets.add(j,\/\*(1.+jetUncert\[j\])\*\/jetCorrL1\[j\]/g' CODE/applyCuts.C

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