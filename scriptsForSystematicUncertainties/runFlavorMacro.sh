#!/bin/bash
start=$(date +%s) ## get current time/date to $start

cd ../

# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the Flavor  uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) Deactivate trigger requirement
for i in CODE/myDeclarations.h ; sed -i 's/applyTriggeronMC=true;/applyTriggeronMC=false;/g' "$i"
# 2) Deactivate PU weighting
for i in CODE/myDeclarations.h ; sed -i 's/PUreweighting=true;/PUreweighting=false;/g' "$i"
# 3) Activate flavor Uncertainty evaluation -Algo
for i in CODE/myDeclarations.h ; sed -i 's/const bool flavorUncertaintyEvaluationAlgo=false;/const bool flavorUncertaintyEvaluationAlgo=true;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const bool flavorUncertaintyEvaluationPhys=true;/const bool flavorUncertaintyEvaluationAlgo=false;/g' "$i"
#
# ----------------------------------------------------------------------------------------------------

# Run whole analysis for MC (together - algo)
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/const int setFlavorSelection=.;/const int setFlavorSelection=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_algo ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_algo plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_algo_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_algo
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (quarks -algo )
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/const int setFlavorSelection=.;/const int setFlavorSelection=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_algo ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_algo plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_algo_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_algo
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (gluons - algo)
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/const int setFlavorSelection=.;/const int setFlavorSelection=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_algo ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_algo plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_algo_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_algo
# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# 1) Activate flavor Uncertainty evaluation - phys
for i in CODE/myDeclarations.h ; sed -i 's/const bool flavorUncertaintyEvaluationAlgo=true;/const bool flavorUncertaintyEvaluationAlgo=false;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const bool flavorUncertaintyEvaluationPhys=false;/const bool flavorUncertaintyEvaluationPhys=true;/g' "$i"
#
# ----------------------------------------------------------------------------------------------------

# Run whole analysis for MC (together - phys)
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/const int setFlavorSelection=.;/const int setFlavorSelection=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_phys ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_phys plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_phys_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/together_phys
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (quarks - phys)
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/const int setFlavorSelection=.;/const int setFlavorSelection=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_phys ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_phys plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_phys_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/quarks_phys
# ----------------------------------------------------------------------------------------------------
# Run whole analysis for MC (gluons - phys)
# ----------------------------------------------------------------------------------------------------
for i in CODE/myDeclarations.h ; sed -i 's/const int setFlavorSelection=.;/const int setFlavorSelection=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_phys ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_phys plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_phys_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_FlavorUncertainty/gluons_phys
# ----------------------------------------------------------------------------------------------------

# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) Activate trigger requirement
for i in CODE/myDeclarations.h ; sed -i 's/applyTriggeronMC=false;/applyTriggeronMC=true;/g' "$i"
# 2) Activate PU weighting
for i in CODE/myDeclarations.h ; sed -i 's/PUreweighting=false;/PUreweighting=true;/g' "$i"
# 3) Deactivate flavor uncertainty evaluation
for i in CODE/myDeclarations.h ; sed -i 's/const bool flavorUncertaintyEvaluationPhys=true;/const bool flavorUncertaintyEvaluationPhys=false;/g' "$i"

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