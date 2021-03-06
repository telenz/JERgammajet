#!/bin/bash
start=$(date +%s) ## get current time/date to $start

cd ../

# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the QCD uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) Switch to ak7PFCHS jets
sed -i 's/jetType=2;/jetType=4;/g' CODE/myDeclarations.h
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
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_mc ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_mc_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_mc plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_mc_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_mc
# ----------------------------------------------------------------------------------------------------

# Run whole analysis for Data
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_data.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_data ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_data_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_data plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_data_SAVED
fi
cp -r plots_2012/PF_L1CHS/data/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak7PFCHS_data
# ----------------------------------------------------------------------------------------------------

# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) Switch to ak7PFCHS jets
sed -i 's/jetType=4;/jetType=2;/g' CODE/myDeclarations.h

# Run whole analysis for MC (ak5)
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_mc ]
  then
    rm -rf plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_mc_SAVED
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_mc plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_mc_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_mc
# ----------------------------------------------------------------------------------------------------

# Run whole analysis for Data (ak5)
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_data.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_data ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_data plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_data_SAVED
fi
cp -r plots_2012/PF_L1CHS/data/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_MCUncertainty/ak5PFCHS_data



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