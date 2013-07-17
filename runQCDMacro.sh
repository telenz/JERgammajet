#!/bin/bash
start=$(date +%s) ## get current time/date to $start


# ----------------------------------------------------------------------------------------------------
# Change relevant stuff for running the QCD uncertainty
# ----------------------------------------------------------------------------------------------------
# 1) deactivate trigger requirement
for i in CODE/myDeclarations.h ; sed -i 's/applyTriggeronMC=true;/applyTriggeronMC=false;/g' "$i"
# 2) deactivate PU weighting
for i in CODE/myDeclarations.h ; sed -i 's/PUreweighting=true;/PUreweighting=false;/g' "$i"
# 3) Activate QCD uncertainty evaluation
for i in CODE/myDeclarations.h ; sed -i 's/QCDUncertaintyEvaluation=false;/QCDUncertaintyEvaluation=true;/g' "$i"
# 4) Change eta Binning
for i in CODE/myDeclarations.h ; sed -i 's/nEtaBins=4;/nEtaBins=2;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const double etaBins\[nEtaBins+1\]={0.0,0.5,1.1,1.7,2.3};/const double etaBins\[nEtaBins+1\]={0.0,1.1,2.3};/g' "$i"
# 5) Change alpha Binning
for i in CODE/myDeclarations.h ; sed -i 's/nAlphaBins=6;/nAlphaBins=4;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const double alphaBins\[nAlphaBins+1\]={0.0,7.5,10.0,12.5,15.0,17.5,20.0};/const double alphaBins\[nAlphaBins+1\]={0.0,8.0,12.0,16.0,20.0};/g' "$i"
# 6) Change pt Binning
for i in CODE/myDeclarations.h ; sed -i 's/nPtBins=12;/nPtBins=5;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const double ptBins\[nPtBins+1\]={22,36,60,88,105,148.5,165,176,200,250,300,400,1000000};/const double ptBins\[nPtBins+1\]={22,100,200,300,400,1000000};/g' "$i"
# 7) Fix N to zero
for i in CODE/myClasses.C ; sed -i 's/\/\/fResolution -> FixParameter(0,0.);/fResolution -> FixParameter(0,0.);/g' "$i"
# 8) Set length to 3
for i in CODE/myClasses.C ; sed -i 's/if(length>=4){/if(length>=3){/g' "$i"


# ----------------------------------------------------------------------------------------------------

# Run whole analysis for MC (woQCD)
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_mc.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_mc.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
cp -r plots_2012/PF_L1CHS/mc/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/mcWoCD

# Run whole analysis for Data (withQCD)
# ----------------------------------------------------------------------------------------------------
root -b -l -q jetphoton_data.C+"(10**9,1)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
root -b -l -q jetphoton_data.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=2;/g' "$i"
root -b -l -q jetphoton_data.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=1;/g' "$i"
root -b -l -q jetphoton_data.C+"(10**9,2)"
for i in CODE/myDeclarations.h ; sed -i 's/detJER=.;/detJER=3;/g' "$i"
cp -r plots_2012/PF_L1CHS/data/root_files plots_2012/PF_L1CHS/systematicUncertainties/scripts/root_files_QCDUncertainty/dataWithQCD
# ----------------------------------------------------------------------------------------------------


# Make changes back
# ----------------------------------------------------------------------------------------------------
# 1) deactivate trigger requirement
for i in CODE/myDeclarations.h ; sed -i 's/applyTriggeronMC=false;/applyTriggeronMC=true;/g' "$i"
# 2) deactivate PU weighting
for i in CODE/myDeclarations.h ; sed -i 's/PUreweighting=false;/PUreweighting=true;/g' "$i"
# 3) Deactivate QCD uncertainty evaluation
for i in CODE/myDeclarations.h ; sed -i 's/QCDUncertaintyEvaluation=true;/QCDUncertaintyEvaluation=false;/g' "$i"
# 4) Change eta Binning
for i in CODE/myDeclarations.h ; sed -i 's/nEtaBins=2;/nEtaBins=4;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const double etaBins\[nEtaBins+1\]={0.0,1.1,2.3};/const double etaBins\[nEtaBins+1\]={0.0,0.5,1.1,1.7,2.3};/g' "$i"
# 5) Change alpha Binning
for i in CODE/myDeclarations.h ; sed -i 's/nAlphaBins=4;/nAlphaBins=6;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const double alphaBins\[nAlphaBins+1\]={0.0,8.0,12.0,16.0,20.0};/const double alphaBins\[nAlphaBins+1\]={0.0,7.5,10.0,12.5,15.0,17.5,20.0};/g' "$i"
# 6) Change pt Binning
for i in CODE/myDeclarations.h ; sed -i 's/nPtBins=5;/nPtBins=12;/g' "$i"
for i in CODE/myDeclarations.h ; sed -i 's/const double ptBins\[nPtBins+1\]={22,100,200,300,400,1000000};/const double ptBins\[nPtBins+1\]={22,36,60,88,105,148.5,165,176,200,250,300,400,1000000};/g' "$i"
# 6) Fix N to zero
for i in CODE/myClasses.C ; sed -i 's/fResolution -> FixParameter(0,0.);/\/\/fResolution -> FixParameter(0,0.);/g' "$i"
# 8) Set length to 4
for i in CODE/myClasses.C ; sed -i 's/if(length>=3){/if(length>=4){/g' "$i"

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