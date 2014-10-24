#!/bin/bash
start=$(date +%s) ## get current time/date to $start

cd ../

# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the QCD uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) deactivate trigger requirement
sed -i 's/applyTriggeronMC=true;/applyTriggeronMC=false;/g' CODE/myDeclarations.h
# 2) deactivate PU weighting
sed -i 's/PUreweighting=true;/PUreweighting=false;/g' CODE/myDeclarations.h
# 3) Activate QCD uncertainty evaluation
sed -i 's/QCDUncertaintyEvaluation=false;/QCDUncertaintyEvaluation=true;/g' CODE/myDeclarations.h
# 4) Change eta Binning
sed -i 's/nEtaBins=4;/nEtaBins=2;/g' CODE/myDeclarations.h
sed -i 's/const double etaBins\[nEtaBins+1\]={0.0,0.5,1.1,1.7,2.3};/const double etaBins\[nEtaBins+1\]={0.0,1.3,2.3};/g' CODE/myDeclarations.h
# 5) Change alpha Binning
sed -i 's/nAlphaBins=6;/nAlphaBins=4;/g' CODE/myDeclarations.h
sed -i 's/const double alphaBins\[nAlphaBins+1\]={0.0,7.5,10.0,12.5,15.0,17.5,20.0};/const double alphaBins\[nAlphaBins+1\]={0.0,8.0,12.0,16.0,20.0};/g' CODE/myDeclarations.h
# 6) Change pt Binning
sed -i 's/nPtBins=12;/nPtBins=5;/g' CODE/myDeclarations.h
sed -i 's/const double ptBins\[nPtBins+1\]={22,36,60,88,105,148.5,165,176,200,250,300,400,1000000};/const double ptBins\[nPtBins+1\]={22,100,200,300,400,1000000};/g' CODE/myDeclarations.h
# 7) Fix N to zero
sed -i 's/\/\/fResolution -> FixParameter(0,0.);/fResolution -> FixParameter(0,0.);/g' CODE/myClasses.C
# 8) Set length to 3
sed -i 's/if(length>=4){/if(length>=3){/g' CODE/myDeclarations.h


# ----------------------------------------------------------------------------------------------------

# Run whole analysis for MC (woQCD)
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_mc.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_mc.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/mcWoQCD ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/mcWoQCD plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/mcWoQCD_SAVED
fi
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/mcWoQCD

# Run whole analysis for Data (withQCD)
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_data.C+"(10**9,1)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=2;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=1;/g' CODE/myDeclarations.h
root -b -l -q jetphoton_data.C+"(10**9,2)"
sed -i 's/detJER=.;/detJER=3;/g' CODE/myDeclarations.h
if [ -d plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/dataWithQCD ]
  then
    mv plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/dataWithQCD plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/dataWithQCD_SAVED
fi
cp -r plots_2012/PF_L1CHS/data/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/dataWithQCD
# ----------------------------------------------------------------------------------------------------


# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) deactivate trigger requirement
sed -i 's/applyTriggeronMC=false;/applyTriggeronMC=true;/g' CODE/myDeclarations.h
# 2) deactivate PU weighting
sed -i 's/PUreweighting=false;/PUreweighting=true;/g' CODE/myDeclarations.h
# 3) Deactivate QCD uncertainty evaluation
sed -i 's/QCDUncertaintyEvaluation=true;/QCDUncertaintyEvaluation=false;/g' CODE/myDeclarations.h
# 4) Change eta Binning
sed -i 's/nEtaBins=2;/nEtaBins=4;/g' CODE/myDeclarations.h
sed -i 's/const double etaBins\[nEtaBins+1\]={0.0,1.3,2.3};/const double etaBins\[nEtaBins+1\]={0.0,0.5,1.1,1.7,2.3};/g' CODE/myDeclarations.h
# 5) Change alpha Binning
sed -i 's/nAlphaBins=4;/nAlphaBins=6;/g' CODE/myDeclarations.h
sed -i 's/const double alphaBins\[nAlphaBins+1\]={0.0,8.0,12.0,16.0,20.0};/const double alphaBins\[nAlphaBins+1\]={0.0,7.5,10.0,12.5,15.0,17.5,20.0};/g' CODE/myDeclarations.h
# 6) Change pt Binning
sed -i 's/nPtBins=5;/nPtBins=12;/g' CODE/myDeclarations.h
sed -i 's/const double ptBins\[nPtBins+1\]={22,100,200,300,400,1000000};/const double ptBins\[nPtBins+1\]={22,36,60,88,105,148.5,165,176,200,250,300,400,1000000};/g' CODE/myDeclarations.h
# 6) Fix N to zero
sed -i 's/fResolution -> FixParameter(0,0.);/\/\/fResolution -> FixParameter(0,0.);/g' CODE/myClasses.C
# 8) Set length to 4
sed -i 's/if(length>=3){/if(length>=4){/g' CODE/myDeclarations.h

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