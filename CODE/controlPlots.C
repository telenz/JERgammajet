#include "myDeclarations.h"
#include "myFunctions.h"
#include "myClasses.h"
#include "TVector2.h"
#include "TMath.h"


TH1D* hPhotonPt[numTrigger]    = {0};
TH1D* hVtxPtbinned[numTrigger] = {0};
TH1D* hVtxPtbinnedwoWeights[numTrigger] = {0};
TH1D* hPUgenMC[numTrigger]           = {0};
TH1D* hPhoton1Pt          = 0;  
TH1D* hPhoton1PtwoWeights = 0;  
TH1D* hJet1Pt       = 0;  
TH1D* hJet2Pt       = 0;
TH1D* hDeltaPhi     = 0;
TH1D* hNPV          = 0;
TH1D* hNPU              = 0;
TH1D* hInvMassMETPhoton = 0;
TH1D* hWeights = 0;

TH1D* hMet             = 0; 
TH2D* hMetPhiPhotonPhi = 0;
TH2D* hWeightsPhotonPt = 0;

// Photon Isolation Variables
TH1D* hEcalRecHit      = 0;
TH1D* hSigmaIetaIeta   = 0;
TH1D* hHcalTowerSum    = 0;
TH1D* hTrkHollowCone   = 0;
TH1D* hHoverE          = 0;
TH1D* pixelSeedVeto    = 0;

TH2D* hPhotonPtJetEtaFlavorGluon_phys       = 0;
TH2D* hPhotonPtJetEtaFlavorLightQuarks_phys = 0;
TH2D* hPhotonPtJetEtaFlavorBottom_phys      = 0;
TH2D* hPhotonPtJetEtaFlavorCharm_phys       = 0;
TH2D* hPhotonPtJetEtaFlavorNonDefined_phys  = 0;

TH2D* hPhotonPtJetEtaFlavorGluon_algo       = 0;
TH2D* hPhotonPtJetEtaFlavorLightQuarks_algo = 0;
TH2D* hPhotonPtJetEtaFlavorBottom_algo      = 0;
TH2D* hPhotonPtJetEtaFlavorCharm_algo       = 0;
TH2D* hPhotonPtJetEtaFlavorNonDefined_algo  = 0;


TH1D* hChiSquareIntrinsic  = 0;
TH1D* hChiSquareImbalance  = 0;
TH1D* hChiSquareResolution = 0;

void initializeControlPlots(){

  char tname[500];
  TString histoName;

  for(int i=0; i<numTrigger; i++){

    if(i<numTrigger-1) sprintf(tname,"Photon Pt for %4.1f GeV < p_{T} < %4.1f GeV ",ptBins[i],ptBins[i+1]);
    else    sprintf(tname,"Photon Pt for %4.1f GeV < p_{T}",ptBins[i]);
    histoName.Form("hPhotonPt%i",i);
    hPhotonPt[i] = createTH1(hPhotonPt[i],histoName,tname,2000,0,2000,"Photon Pt");

    if(i<nPtBins-1) sprintf(tname,"#Vtx %4.1f GeV < p_{T} < %4.1f GeV ",ptBins[i],ptBins[i+1]);
    else           sprintf(tname,"#Vtx for %4.1f GeV < p_{T}",ptBins[i]);
    histoName.Form("hVtxPtBinned%i",i);
    hVtxPtbinned[i] = createTH1(hVtxPtbinned[i],histoName,tname,60,0,60,"Number of Vertices");
    histoName.Form("hVtxPtBinnedwoWeights%i",i);
    hVtxPtbinnedwoWeights[i] = createTH1(hVtxPtbinnedwoWeights[i],histoName,tname,60,0,60,"Number of Vertices");
  
    if(i<numTrigger-1) sprintf(tname,"PU distribution MC for %4.1f GeV < p_{T} < %4.1f GeV ",ptBins[i],ptBins[i+1]);
    else               sprintf(tname,"PU distribution MC for %4.1f GeV < p_{T}",ptBins[i]);
    histoName.Form("hPUgenMC%i",i);
    hPUgenMC[i] = createTH1(hPUgenMC[i],histoName,tname,600,0,60,"#PU");  
  }

  
  hPhoton1Pt          = createTH1(hPhoton1Pt,"hPhoton1Pt","leading Photon p_{T} (after all cuts)",3000,0,3000,"p_{T}^{#gamma} (GeV)");
  hPhoton1PtwoWeights = createTH1(hPhoton1PtwoWeights,"hPhoton1PtwoWeights","leading Photon p_{T} (after all cuts)",3000,0,3000,"p_{T}^{#gamma} (GeV)");
  hJet1Pt    = createTH1(hJet1Pt,"hJet1Pt","Jet12 p_{T}",2000,0,2000,"p_{T} (GeV)");
  hJet2Pt    = createTH1(hJet2Pt,"hJet2Pt","Jet p_{T}",2000,0,2000,"p_{T} (GeV)");
  hDeltaPhi  = createTH1(hDeltaPhi,"hDeltaPhi","#Delta #Phi",100,0,3.142,"#Delta #Phi");
  hNPV       = createTH1(hNPV,"hNPV","Number of Vertices",60,0,60,"Number of vertices");
  hNPU       = createTH1(hNPU,"hNPU","NPU",60,0,60,"NPU");
  hMet       = createTH1(hMet,"hMet","Met",1000,0,1000,"Met");
  hInvMassMETPhoton = createTH1(hInvMassMETPhoton,"hInvMassMETPhoton","hInvMassMETPhoton",1000,0,1000,"inv. mass");
  double xBins[7] = {0, 0.0003, 0.002,0.01, 0.4,23,83};
  hWeights         = new TH1D("hWeights","hWeights",6,xBins);
  hWeightsPhotonPt = new TH2D("hWeightsPhotonPt","hWeightsPhotonPt",6,xBins,nPtBins,ptBins);

  hMetPhiPhotonPhi = createTH2(hMetPhiPhotonPhi,"hMetPhiPhotonPhi",100,0,3.2,100,0,3.2,"Photon_Phi","MET_Phi");

  hEcalRecHit     = createTH1(hEcalRecHit,"hEcalRecHit","hEcalRecHit",1000,0,19,"EcalRecHit");
  hSigmaIetaIeta  = createTH1(hSigmaIetaIeta,"hSigmaIetaIeta","hSigmaIetaIeta",1000,0,0.04,"SigmaIetaIeta");
  hHcalTowerSum   = createTH1(hHcalTowerSum,"hHcalTowerSum","hHcalTowerSum",1000,0,13,"HcalTowerSum");
  hTrkHollowCone  = createTH1(hTrkHollowCone,"hTrkHollowCone","hTrkHollowCone",1000,0,4,"TrkHollowCone");
  hHoverE         = createTH1(hHoverE,"hHoverE","hHoverE",1000,0,0.6,"HoverE");
  pixelSeedVeto   = createTH1(pixelSeedVeto,"pixelSeedVeto","pixelSeedVeto",2,0,2,"pixelSeed");
  hPhotonPtJetEtaFlavorGluon_phys       = createTH2(hPhotonPtJetEtaFlavorGluon_phys,"hPhotonPtJetEtaFlavorGluon_phys",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorLightQuarks_phys = createTH2(hPhotonPtJetEtaFlavorLightQuarks_phys,"hPhotonPtJetEtaFlavorLightQuarks_phys",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorBottom_phys      = createTH2(hPhotonPtJetEtaFlavorBottom_phys,"hPhotonPtJetEtaFlavorBottom_phys",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorCharm_phys       = createTH2(hPhotonPtJetEtaFlavorCharm_phys,"hPhotonPtJetEtaFlavorCharm_phys",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorNonDefined_phys  = createTH2(hPhotonPtJetEtaFlavorNonDefined_phys,"hPhotonPtJetEtaFlavorNonDefined_phys",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorGluon_algo       = createTH2(hPhotonPtJetEtaFlavorGluon_algo,"hPhotonPtJetEtaFlavorGluon_algo",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorLightQuarks_algo = createTH2(hPhotonPtJetEtaFlavorLightQuarks_algo,"hPhotonPtJetEtaFlavorLightQuarks_algo",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorBottom_algo      = createTH2(hPhotonPtJetEtaFlavorBottom_algo,"hPhotonPtJetEtaFlavorBottom_algo",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorCharm_algo       = createTH2(hPhotonPtJetEtaFlavorCharm_algo,"hPhotonPtJetEtaFlavorCharm_algo",1000,0,1000,500,0,5,"photonPt","jetEta");
  hPhotonPtJetEtaFlavorNonDefined_algo  = createTH2(hPhotonPtJetEtaFlavorNonDefined_algo,"hPhotonPtJetEtaFlavorNonDefined_algo",1000,0,1000,500,0,5,"photonPt","jetEta");

  
}


void initializeControlPlotsInCalcScale(){


  hChiSquareIntrinsic  = createTH1(hChiSquareIntrinsic,"hChiSquareIntrinsic","hChiSquareIntrinsic",100,0,30,"#Chi^{2}/ndof");
  hChiSquareImbalance  = createTH1(hChiSquareImbalance,"hChiSquareImbalance","hChiSquareImbalance",100,0,30,"#Chi^{2}/ndof");
  hChiSquareResolution = createTH1(hChiSquareResolution,"hChiSquareResolution","hChiSquareResolution",100,0,30,"#Chi^{2}/ndof");
  
}

void plotControlPlots(){


  cout<<endl<<"Plot Control Plots!!"<<endl<<endl;

  TString texte;

  for(int i=0; i<numTrigger; i++){
    texte.Form("VtxN%i",i);
    plotTH1(hVtxPtbinned[i],texte,0); 

    texte.Form("VtxNwoWeights%i",i);
    plotTH1(hVtxPtbinnedwoWeights[i],texte,0);   

    texte.Form("PUgenMC%i",i);
    plotTH1(hPUgenMC[i],texte,0);
  }

  plotTH1(hPhoton1Pt,"hPhoton1Pt",0);
  plotTH1(hPhoton1PtwoWeights,"hPhoton1PtwoWeights",0);
  plotTH1(hJet1Pt,"hJet1Pt",0);
  plotTH1(hJet2Pt,"hJet2Pt",0);
  plotTH1(hDeltaPhi,"hDeltaPhi",0);
  plotTH1(hMet,"hMet",0);
  plotTH2(hMetPhiPhotonPhi,"hMetPhiPhotonPhi");
  plotTH2(hWeightsPhotonPt,"hWeightsPhotonPt");
  plotTH1(hWeights,"hWeights",0);

  //plotTH1(,"",0);
  plotTH1(hEcalRecHit,"hEcalRecHit",0);
  plotTH1(hSigmaIetaIeta,"hSigmaIetaIeta",0);
  plotTH1(hHcalTowerSum,"hHcalTowerSum",0);
  plotTH1(hTrkHollowCone,"hTrkHollowCone",0);
  plotTH1(hHoverE,"hHoverE",0);
  plotTH1(pixelSeedVeto,"hpixelSeedVeto",0);

  if(isMC){
    plotTH2(hPhotonPtJetEtaFlavorGluon_phys,"hPhotonPtJetEtaFlavorGluon_phys");
    plotTH2(hPhotonPtJetEtaFlavorLightQuarks_phys,"hPhotonPtJetEtaFlavorLightQuarks_phys");
    plotTH2(hPhotonPtJetEtaFlavorBottom_phys,"hPhotonPtJetEtaFlavorBottom_phys");
    plotTH2(hPhotonPtJetEtaFlavorCharm_phys,"hPhotonPtJetEtaFlavorCharm_phys");
    plotTH2(hPhotonPtJetEtaFlavorNonDefined_phys,"hPhotonPtJetEtaFlavorNonDefined_phys"); 
    plotTH2(hPhotonPtJetEtaFlavorGluon_algo,"hPhotonPtJetEtaFlavorGluon_algo");
    plotTH2(hPhotonPtJetEtaFlavorLightQuarks_algo,"hPhotonPtJetEtaFlavorLightQuarks_algo");
    plotTH2(hPhotonPtJetEtaFlavorBottom_algo,"hPhotonPtJetEtaFlavorBottom_algo");
    plotTH2(hPhotonPtJetEtaFlavorCharm_algo,"hPhotonPtJetEtaFlavorCharm_algo");
    plotTH2(hPhotonPtJetEtaFlavorNonDefined_algo,"hPhotonPtJetEtaFlavorNonDefined_algo"); 
  }
  if(isMC){
    plotTH1(hNPU,"hNPU",0);
    plotTH1(hNPV,"hNPV",0);
  }
}


void plotControlPlotsInCalcScale(){


  cout<<endl<<"Plot Control Plots!!"<<endl<<endl;


  if(isMC){
    plotTH1(hChiSquareIntrinsic,"hChiSquareIntrinsic",0);
    plotTH1(hChiSquareImbalance,"hChiSquareImbalance",0);
  }
  plotTH1(hChiSquareResolution,"hChiSquareResolution",0);
 

}
