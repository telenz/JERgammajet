#ifndef myDeclarations_h
#define myDeclarations_h

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TString.h"
#include "TLorentzVector.h"

// Numbers of alpha, pTGamma and eta bins 
const int nPtBins=12;      // numbers of Pt bins
const int nAlphaBins=6;    // numbers of alpha bins 
const int nEtaBins=4;      // number of eta bins

// Definition of all bin intervals
const double ptBins[nPtBins+1]={22,36,60,88,105,148.5,165,176,200,250,300,400,1000000};
const double alphaBins[nAlphaBins+1]={0.0,7.5,10.0,12.5,15.0,17.5,20.0};
const double etaBins[nEtaBins+1]={0.0,0.5,1.1,1.7,2.3};


const int numTrigger = 7;   // 8 for 2012 and 5 for 2011

const int detJER=3;         // = 1 (2sigma Gaus) = 2 (RMS 95%) = 3 (RMS 99%)
const int jetType=2;        // 1 = PF_L1Fast , 2 = PF_CHS Jets , 3 = Calo_L1FastJ, 4 = ak7 PF_CHS Jetset 
const int date   = 2012;    // 2012 or 2011
const int set    = 5;       // 1=AB; 2= ABC; 3=AB rereco; 4=C; 5=ABCD rereco

const bool applyTriggeronMC=true;  // shall the triggers also applied on MC
const bool PUreweighting=true;     // shall MC be PU reweighted to the data distribution
const bool addQCD=false;           // add QCD sample
const bool addWJet=false;          // add W+jet sample

const bool testClosure=false;                         // for ratio closure test
const bool QCDUncertaintyEvaluation=false;             // Set on true if QCD sample shall be added to the data procedure
const bool flavorUncertaintyEvaluationPhys=false;     //Set on true if you want to evaluate flavor Uncertainty 
const bool flavorUncertaintyEvaluationAlgo=false;     //Set on true if you want to evaluate flavor Uncertainty 
const int setFlavorSelection=3;                       //1 = together, 2 = all quarks, 3 = gluons, 4 = uds, 5 = cb 

const double deltaPhiCutValue=2.95;  // Cut value of deltaPhi

// Declaration of the path where all files shall be saved in the end 
TString PDFPath, RootPath, DataType, MethodType; 
     
// To save strange events for 2011
ofstream filestr, EventVariables, cutflow;

TChain* chain;

// sum indices for the jet collection
int idx2ndJet = -1;
int idx1stJet = -1;
int idxPhoton =-1;
int gen1stJetidx = -1;
int gen2ndJetidx = -1;

// Define some global variables 
float response  = 0;
float imbalance = 0;
float intrinsic = 0; 
float alpha     = 0; 

// 1st and 2nd Jet Pt
float jetPt2nd    = 0;
float jetPt1stJet = 0;

TLorentzVector photonVector;
TLorentzVector jet2ndVector;
TLorentzVector jet1stVector;
TLorentzVector sumVector;

float PUWeight = 1;
int count1stJetNotMatched = 0;
int no2ndJetinEvent = 0;


// ------ Book variables to be read from tree [needed in: readGammaJet.C] ------------------
const int NJETS    = 200;   // Maximum number of jets per event
const int NPHOTONS = 100;   // Maximum number of Photons per event

float crossSection;
int nobjPhoton;                // Number of Photons in an event
int nobjJet;                   // Number of jets in an event
bool tight[NPHOTONS];          // tight Photon or not
bool loose[NPHOTONS];          // loose Photon or not
bool jetIDTight[NJETS];         // tight jet or not
bool jetIDLoose[NJETS];         // loose jet or not
float photonPt[NPHOTONS];      // Pt of Photons in an event
float photonE[NPHOTONS];      // E of Photons in an event
float photonEta[NPHOTONS];     // Photon Eta
float photonPhi[NPHOTONS];     // Phi of Photons in an event  
float genPhotonPt[NPHOTONS];
float genPhotonEta[NPHOTONS];     // Photon Eta
float genPhotonPhi[NPHOTONS];     // Phi of Photons in an event 
float jetPt[NJETS];            // Pt of jets in an event
float jetCorrL2L3[NJETS];	 // Jet energy correction factors (L2 = flat vs eta  ;L3 = flat vs pt; for data also residual corrections are included)
float jetCorrL1[NJETS];        // Jet energy correction Factor L1 (pileup corrections)
float jetUncert[NJETS];        // Jet energy correction uncertainty
float jetPhi[NJETS];           // Phi of jets in an event
float jetEta[NJETS];	         // Jet eta  
float jetE[NJETS];          // Generated Jet Pt

int vtxN;                     // number of vertices in an event 
bool EcalBEFilter; //Ecal bool BE Filter = Boundary Bnergy Filter
bool EcalTPFilter; //Ecal bool TP Filter = Trigger Primitive Filter
unsigned int runNum;    // run number
unsigned int lumiNum;   // lumi block number 
unsigned int eventNum;  // event number 
float photonIsoEcal[NPHOTONS];   // isolation Photon Ecal
float photonIsoHcal[NPHOTONS];   // isolation Photon Hcal
float photonIsoTrk[NPHOTONS];    // isolation Photon Tracker 
float photonSigmaIetaIeta[NPHOTONS];    
float photonHadronicOverEM[NPHOTONS];
float jetFRBX[NJETS];            // fraction of jet energy carried by the most energetic readout box
float jetFHPD[NJETS];            // fraction of jet energy carried by the ``hottest'' (or most energetic) HPD (Hybrid Photo Diodes)
float jetEMF[NJETS];             // energy fraction measured in the calorimeter (E_ECAL/E_{HCAL+ECAL})/ EMF = ElectroMagnetic Fraction
float weight;
// only needed in in calcPUWeights.C
float PUMCNumTruth;            // Number of PU events
// MC specific
float genJetPt[NJETS];          // Generated Jet Pt
float genJetEta[NJETS];          // Generated Jet Pt
float genJetPhi[NJETS];          // Generated Jet Pt
bool photonHasPixelSeed[NPHOTONS];
float genEvtScale;

float met_T1;
float metPhi_T1;

int genJetID_phys[NJETS];            // Generated Jet ID
int genJetID_algo[NJETS];            // Generated Jet ID
float rho;
bool hltPhoton[numTrigger];     // Trigger
float photonEt;
unsigned int genJetColJetIdx[NJETS];

int nobjGenJet;
int nobjGenPhoton;

int nMax;
//--------------------------------------------------------------------------------------------------------------------------------------------
// Cut variables for the cutflow
int cut0=0;
int cut1=0;
int cut2=0;
int cut3=0;
int cut4=0;
int cut5=0;
int cut6=0;
int cut7=0;
int cut8=0;
int cut9=0;
int cut10=0;
int cut11=0; 
int cut12[numTrigger] = {0};
int cut13[3] = {0};
int cut14 = 0;
int cut15 = 0;
int nocut=0;  

//------------------------ Histograms ----------------------------------------
float deltaRPhoton1stJet = 0;
float deltaRPhoton2ndJet = 0;
float deltaR1stJet2ndJet = 0;

TH1D* hDeltaPhi1st2ndJet = 0;
TH2D* hDeltaPhi1st2ndJetDeltaPt = 0;
TH2D* hJet1PtPhoton1Pt   = 0;
TH2D* hImbalanceJetN     = 0;   
TH2D* hImbalancePhoton1Pt= 0;    
TH2D* hImbalanceVtx      = 0;    
TH2D* hImbalanceJetPt    = 0;     
TH2D* hImbalanceJetEta   = 0;    
TH2D* hImbalancePhotonEta=0;
TH2D* hImbalanceWeights = 0;
TH2D* hImbalanceJet2Pt =0;
TH2D* hImbalanceGenJet2Pt =0;
TH2D* hJet2PtGenJet2Pt =0;
TH2D* hImbalanceAlpha = 0;
TH2D* hImbalanceDeltaRPhoton1stJet =0;
TH2D* hImbalanceDeltaRPhoton2ndJet =0;
TH2D* hImbalanceDeltaR1stJet2ndJet =0;
TH2D* hImbalanceEta2ndJet = 0;
TH2D* hImbalanceRatioPhotonPt12 = 0;
TH2D* hImbalanceRatioPhotongenPhoton = 0;
TH2D* hImbalanceDeltaPhiPhoton2ndJet = 0;
TH2D* hImbalancePhoton2Pt = 0;
TH2D* hImbalanceDeltaEtaPhoton2ndJet = 0;
TH2D* hImbalanceNobjPhoton = 0;
TH2D* hImbalance2ndJettight = 0;
TH2D* hImbalance2ndJetCorr = 0;
TH2D* hImbalance1stGenJetID = 0;
TH2D* hImbalance2ndGenJetID = 0;
TH2D* hPhotonEta2ndJet = 0;
TH1D* hEnergyBalance = 0;
TH2D* hImbalanceRatioPhotonPtJetColPt = 0;
TH1D* h2ndJetPt1stJetHemisphere[nAlphaBins] = {0};
TH1D* h2ndJetPtPhotonHemisphere[nAlphaBins] = {0};
TH1D* h2ndJetEta1stJetHemisphere[nAlphaBins] = {0};
TH1D* h2ndJetEtaPhotonHemisphere[nAlphaBins] = {0};
TH1D* hAlpha = 0;
TH2D* hJet1PtPhoton1Pt2 =0;
TH2D* hImbalanceJetN2     = 0;   
TH2D* hImbalancePhoton1Pt2= 0;    
TH2D* hImbalanceVtx2   = 0;    
TH2D* hImbalanceJetPt2   = 0;     
TH2D* hImbalanceJetEta2   = 0;    
TH2D* hImbalancePhotonEta2 =0;
TH2D* hImbalanceWeights2 = 0;
TH2D* hImbalanceJet2Pt2 = 0;
TH2D* hImbalanceGenJet2Pt2 =0;
TH2D* hJet2PtGenJet2Pt2 =0;
TH2D* hImbalanceAlpha2 = 0;
TH2D* hImbalanceDeltaRPhoton1stJet2 = 0;
TH2D* hImbalanceDeltaRPhoton2ndJet2 = 0;
TH2D* hImbalanceDeltaR1stJet2ndJet2 = 0;
TH2D* hImbalanceEta2ndJet2 = 0;
TH2D* hImbalanceRatioPhotonPt122 = 0;
TH2D* hImbalanceRatioPhotongenPhoton2 = 0;
TH2D* hImbalanceDeltaPhiPhoton2ndJet2 = 0;
TH2D* hImbalancePhoton2Pt2 = 0;
TH2D* hImbalanceDeltaEtaPhoton2ndJet2 = 0;
TH2D* hImbalanceNobjPhoton2 = 0;
TH2D* hImbalance2ndJettight2 = 0;
TH2D* hImbalance2ndJetCorr2 = 0;
TH2D* hImbalance1stGenJetID2 = 0;
TH2D* hImbalance2ndGenJetID2 = 0;
TH2D* hPhotonEta2ndJet2 = 0;
TH1D* hJetN         = 0;     // Number of Jets in one event
TH1D* hJetEta       = 0;     // Eta of first Jet in an event
TH1D* hJetEta_1     = 0;     // Eta of first Jet in an event with vtxN>1
TH2D* hJetEtaPt1    = 0;     // 2d eta against Pt of first Jet
TH1D* hEnergyBalance2 = 0;
TH2D* hImbalanceRatioPhotonPtJetColPt2 = 0;
TH1D* PhotonPtwithWeight = 0;
TH1D* PhotonPtwithoutWeight = 0;

TH2D* hResponseLeadJet   = 0;     
TH1D* hPhotonEta_low_resp_reg = 0;
TH2D* hSmallJetPtResponse = 0;
TH2D* hPhotonEta_high_pt_reg_Response = 0;
TH2D* hResponseEnergyFraction=0;
TH2D* hResponseVecEnergyFraction=0;
TH2D* hResponseL1=0;
TH2D* hResponseL2L3=0;
TH2D* hResponseJetEMF = 0;
TH2D* hResponseJetFHPD = 0;
TH2D* hResponseJetFRBX =0; 
TH2D* hResponsePhotonEMF = 0;
TH2D* hResponsePhotonFHPD = 0;
TH2D* hResponsePhotonFRBX =0; 
TH2D* hResponsePhotonEta =0;
TH2D* hResponsePhotonPtRatio = 0;
TH1D* hWeight = 0;
TH1D* hWeightWeight = 0;
TH1D* hPhotonT1Pt   = 0;     // Pt of tight first Photon in two different eta regions
TH1D* hPhotonN      = 0;     // Number of Photons in one event
TH1D* hPhotonNtight = 0;     // Number of tight Photons in one event
TH1D* hPhotonPhi    = 0;     // Phi of first Photon in an event
TH1D* hPhotonEta    = 0;     // Eta of first Photon in an event
TH2D* hPhoton12Pt   = 0;     // 2d plot of Photon pt
TH1D* hVtxN    = 0;
TH1D* hVtxN1[nPtBins][nEtaBins]    = {{0}};
TH1D* hVtxN2[nPtBins][nEtaBins]    = {{0}};
TH1D* hVtxN3[nPtBins][nEtaBins]    = {{0}};
TH1D* hVtxN4[nPtBins][nEtaBins]    = {{0}};
TH1D* hVtxN5[nPtBins][nEtaBins]    = {{0}};
TH1D* hDphi      = 0;        // delta phi of Jet and Photon
TH1D* hRatioPt   = 0;        // Jet2Pt/Photon1Pt
TH2D* hResponsePt= 0;        // 2d Response and Pt of first Photon
TH2D* hResponsePt_0_5_Deta;
TH2D* hResponsePt_1_1_Deta;
TH2D* hResponsePt_1_5_Deta;
TH2D* hResponsePtBeforeEcal= 0;  // 2d Response and Pt of first Photon before ecal cut
TH1D* hPUDataDist[numTrigger] = {0};     //Data Pilup Ditributions
TH1D* hPUMCDist      = 0;     //MC Pileup Ditributions 
TH1D* hPUWeight[numTrigger]   = {0};     //Pileup Weights
TH2D* hPhotonPtVtx   = 0;
TH2D* hPhotonPtVtxBinned[numTrigger-1]   = {0};   // 2d Vertex against Photon Pt for every Trigger withh different prescales
TH2D* hPhotonPtVtxbeforeCUTS = 0;
TH1D* hRho[numTrigger] = {0};
TH1D* hTriggerEffBefore[numTrigger] = {0};
TH1D* hTriggerEffAfter[numTrigger] = {0};
TH1D* hRhoVtxBinned[40] = {0};
TH1D* hPUsysY[4][nPtBins][nEtaBins] = {{{0}}};

#endif
