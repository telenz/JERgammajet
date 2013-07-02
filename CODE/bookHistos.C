#include "myDeclarations.h"
#include "myFunctions.h"

#include "TStyle.h"

void bookHistos() {
 
//-------------------------------- 1d - histograms-------------------------------

  

  // Setting different upper bounds on Pt depending on the tree which is used
  double upp_pt = 1200.;    
  char tname[100]; 

  TString histoName; 


  // what shall be shown in the histogram statistics box
  gStyle->SetOptStat("emr");



  // 1d: Histogram for  Alpha after whole selection
  hAlpha = createTH1(hAlpha,"hAlpha","#alpha",50,0,1,"#alpha");

  

  


  hDeltaPhi1st2ndJet = createTH1(hDeltaPhi1st2ndJet,"hDeltaPhi1st2ndJet","deltaPhi1st2ndJet",100,0,3.2,"deltaPhi1st2ndJet");
  hDeltaPhi1st2ndJetDeltaPt = createTH2(hDeltaPhi1st2ndJetDeltaPt,"hDeltaPhi1st2ndJetDeltaPt",100,0,3.1,100,0,500,"phi","deltaPt");

 
  // 1d: Histogram for jet PT in small response region
  hPhotonEta_low_resp_reg = createTH1(hPhotonEta_low_resp_reg,"hPhotonEtaLowResponseReg","Photon Eta in low response region (response <0.3)",48,-3,3,"#eta^{#gamma}");
 
  // 1d: Histograms for tight Photon Pt 
  hPhotonT1Pt = createTH1(hPhotonT1Pt,"hPhotonT1Pt","Tight Photon1 p_{T}",1200,0,upp_pt,"p_{T,1} (GeV)");     
 
  
  // 1d: Histogram for number of Photons in one Event
  hPhotonN = createTH1(hPhotonN,"hPhotonN","Number of Photons",7,0,7,"# Photons");  

  // 1d: Histogram for number of tight Photons in one Event
  hPhotonNtight = createTH1(hPhotonNtight,"hPhotonNtight","Number of tight Photons",7,0,7,"# tight Photons");

  // 1d: Histogram for number of Jets in one Event
  hJetN = createTH1(hJetN,"hJetN","Number of Jets",60,0,60,"# Jets");

  // 1d: Histogram for Phi of first Photon in an event 
  hPhotonPhi = createTH1(hPhotonPhi,"hPhotonPhi","Phi of first Photon",100,0,3.2,"Phi of first Photon");

  // 1d: Histogram for Eta of first Photon in an event 
  hPhotonEta = createTH1(hPhotonEta,"hPhotonEta","#eta of first Photon",100,-5,5,"#eta_{1}^{#gamma}");
  
  // 1d: Histogram for Eta of first Jet in an event 
  hJetEta = createTH1(hJetEta,"hJetEta","Eta of first Jet",100,-5.5,5.5,"#eta^{Jet_{1}}");
  
  // 1d: Histogram for Eta of first Jet in an event with vtxN>1
  hJetEta_1 = createTH1(hJetEta_1,"hJetEta1","#eta of first Jet with vtx>1",100,-5.5,5.5,"#eta^{Jet_{1}}");
  
  
  hVtxN = createTH1(hVtxN,"hVtxN","Number of Vertices",60,0,60,"Number of vertices");

  for(int i=0; i<nPtBins; i++){
      for(int j=0; j<nEtaBins; j++){
      
	hVtxN1[i][j] = createTH1(hVtxN1[i][j],"","Number of Vertices",60,0,60,"Number of vertices");
	hVtxN2[i][j] = createTH1(hVtxN2[i][j],"","Number of Vertices",60,0,60,"Number of vertices");
	hVtxN3[i][j] = createTH1(hVtxN3[i][j],"","Number of Vertices",60,0,60,"Number of vertices");
	hVtxN4[i][j] = createTH1(hVtxN4[i][j],"","Number of Vertices",60,0,60,"Number of vertices");
	hVtxN5[i][j] = createTH1(hVtxN5[i][j],"","Number of Vertices",60,0,60,"Number of vertices");
      }
  }

  


  
  for(int i=0; i<numTrigger; i++){
    if(i<nPtBins-1) sprintf(tname,"rho %4.1f GeV < p_{T} < %4.1f GeV ",ptBins[i],ptBins[i+1]);
    else    sprintf(tname,"rho for %4.1f GeV < p_{T}",ptBins[i]);
    histoName.Form("hRhoPtBinned%i",i);
    hRho[i] = createTH1(hRho[i],histoName,tname,60,0,60,"#Vtx");
  }
  for(int i=0; i<40; i++){
    sprintf(tname,"rho for #Vtx = %i",i);
    histoName.Form("hRhoVtxBinned%i",i);
    hRhoVtxBinned[i] = createTH1(hRhoVtxBinned[i],histoName,tname,40,0,40,"#Vtx");
  }

  for(int i=0; i<numTrigger; i++){
    if(i<nPtBins-1) sprintf(tname,"#Vtx %4.1f GeV < p_{T} < %4.1f GeV (before)",ptBins[i],ptBins[i+1]);
    else    sprintf(tname,"#Vtx for %4.1f GeV < p_{T} (before)",ptBins[i]);
    histoName.Form("hTriggerEffPtBinnedBefore%i",i);
    hTriggerEffBefore[i] = createTH1(hTriggerEffBefore[i],histoName,tname,60,0,60,"#Vtx");
  }

  for(int i=0; i<numTrigger; i++){
    if(i<nPtBins-1) sprintf(tname,"#Vtx %4.1f GeV < p_{T} < %4.1f GeV (after)",ptBins[i],ptBins[i+1]);
    else    sprintf(tname,"#Vtx for %4.1f GeV < p_{T} (after)",ptBins[i]);
    histoName.Form("hTriggerEffPtBinnedAfter%i",i);
    hTriggerEffAfter[i] = createTH1(hTriggerEffAfter[i],histoName,tname,60,0,60,"#Vtx");
  }
  
 
  // 1d: Histogram for delta Phi between first Jet and first tight Photon 
  hDphi = createTH1(hDphi,"hdphi","#Delta #Phi between first Jet and first tight Photon",100,0.,3.2,"#Delta  #Phi");

  // 1d: Histogram for Pt_Jet2 Pt_Photon Ratio
  hRatioPt = createTH1(hRatioPt,"hRatioPt","PtJet2/PtPhoton1",100,0,3,"Ratio");
 
  //1d weights
  hWeight = createTH1(hWeight,"hWeight","Weights",50,0,0.25,"Weights");
  hWeightWeight = createTH1(hWeightWeight,"hWeightWeight","Weights weighted",50,0,0.25,"Weights"); 

  // 1d PUWeights
  
  for(int i=0; i<numTrigger; i++){
    if(i<numTrigger-1) sprintf(tname,"PU-Weights for %4.1f GeV < p_{T} < %4.1f GeV ",ptBins[i],ptBins[i+1]);
    else    sprintf(tname,"PU-Weights for %4.1f GeV < p_{T}",ptBins[i]);
    histoName.Form("PUWeight%i",i);
    hPUWeight[i] = createTH1(hPUWeight[i],histoName,tname,600,0,60,"# Pileup");
    //hPUWeight[i] = createTH1(hPUWeight[i],histoName,tname,60,0,60,"# Pileup");
  }

  // 
  


  PhotonPtwithWeight = createTH1(PhotonPtwithWeight,"PhotonPtwithWeight","Leading Photon Pt with weight",100,0,2000,"Photon Pt");
  PhotonPtwithoutWeight = createTH1(PhotonPtwithoutWeight,"PhotonPtwithoutWeight","Leading Photon Pt without weight",100,0,2000,"Photon Pt");

  

  // PU Distribution for several PU bins (for sys uncertainty)
  for(int k=0; k<4; k++){
    for(int i=0; i<nPtBins; i++){
      for(int j=0; j<nEtaBins; j++){
	sprintf(tname,"PU Distribution for %d. PU %d Pt and %d eta bin",k+1,i+1,j+1);
	histoName.Form("hPUsysY%iPU%iPT%iEta",k,i,j);
	hPUsysY[k][i][j] = createTH1(hPUsysY[k][i][j],histoName,tname,40,0,40,"PU");
      }
    }
  }
  
  hEnergyBalance =  createTH1(hEnergyBalance,"hEnergyBalance","hEnergyBalance",100,0,2,"sumJetPt/photonPt");
  hEnergyBalance2 =  createTH1(hEnergyBalance2,"hEnergyBalance2","hEnergyBalance2",100,0,2,"sumJetPt/photonPt");
  for(int k=0; k<nAlphaBins; k++){

    histoName.Form("h2ndJetPt1stJetHemisphere_%i_alphaBin",k);
    h2ndJetPt1stJetHemisphere[k] =  createTH1(h2ndJetPt1stJetHemisphere[k],histoName,histoName,400,0,400,"2ndJetPt");
    histoName.Form("h2ndJetPtPhotonHemisphere_%i_alphaBin",k);
    h2ndJetPtPhotonHemisphere[k] =  createTH1(h2ndJetPtPhotonHemisphere[k],histoName,histoName,400,0,400,"2ndJetPt");
    histoName.Form("h2ndJetEta1stJetHemisphere_%i_alphaBin",k);
    h2ndJetEta1stJetHemisphere[k] =  createTH1(h2ndJetEta1stJetHemisphere[k],histoName,histoName,100,-6,6,"2ndJetEta");
    histoName.Form("h2ndJetEtaPhotonHemisphere_%i_alphaBin",k);
    h2ndJetEtaPhotonHemisphere[k] =  createTH1(h2ndJetEtaPhotonHemisphere[k],histoName,histoName,100,-6,6,"2ndJetEta");
  }
  //-------------------------------- 2d - histograms-------------------------------

  
  hJet1PtPhoton1Pt = createTH2(hJet1PtPhoton1Pt,"Jet1 Pt and Photon1 Pt",100,0,upp_pt,100,130,upp_pt,"p_{T}^{Jet}","p_{T}^{#gamma}");
  hImbalanceJetN = createTH2(hImbalanceJetN,"Imbalance against Jet Multiplicity",40,0.4,1.5,100,0,100,"Imbalance","# Jets");
  hImbalancePhoton1Pt = createTH2(hImbalancePhoton1Pt,"Imbalance against Photon1 Pt",40,0.4,1.5,50,0,4000,"Imbalance","Photon1 Pt");
  hImbalanceVtx = createTH2(hImbalanceVtx,"Imbalance against number of vertices",40,0.4,1.5,50,0,50,"Imbalance","# Vtx");
  hImbalanceJetPt= createTH2(hImbalanceJetPt,"Imbalance against jetPt",100,0.4,1.5,100,1,3000,"Imbalance","jetPt");
  hImbalanceJetEta= createTH2(hImbalanceJetEta,"Imbalance against jet eta",100,0.4,1.5,100,0,1.0,"Imbalance","jet Eta");
  hImbalancePhotonEta = createTH2(hImbalancePhotonEta,"Imbalance against #eta_{#gamma}",20,0.4,1.5,100,-1.3,1.3,"Imbalance","#eta_{#gamma}");
  hImbalanceWeights = createTH2(hImbalanceWeights,"Imbalance against Weights",100,0.4,1.5,200,0,TMath::Power(10,-5),"Imbalance","Weights");
  hImbalanceJet2Pt = createTH2(hImbalanceJet2Pt,"Imbalance against Jet2Pt",20,0.4,1.5,500,0,500,"Imbalance","Jet2Pt");
  hImbalanceGenJet2Pt = createTH2(hImbalanceGenJet2Pt,"Imbalance against genJet2Pt",20,0.4,1.5,500,0,500,"Imbalance","genJet2Pt");
  hJet2PtGenJet2Pt = createTH2(hJet2PtGenJet2Pt,"Jet2Pt against genJet2Pt",500,0,500,500,0,500,"Jet2Pt","genJet2Pt");
  hImbalanceAlpha = createTH2(hImbalanceAlpha,"Imbalance against Alpha",80,0.0,2.0,80,0,20,"Imbalance","Alpha");
  hImbalanceDeltaRPhoton1stJet = createTH2(hImbalanceDeltaRPhoton1stJet,"Imbalance against deltaRPhoton1stJet",40,0.4,1.5,50,0,10,"Imbalance","DeltaR");
  hImbalanceDeltaRPhoton2ndJet = createTH2(hImbalanceDeltaRPhoton2ndJet,"Imbalance against deltaRPhoton2ndJet",40,0.4,1.5,50,0,10,"Imbalance","DeltaR");
  hImbalanceDeltaR1stJet2ndJet = createTH2(hImbalanceDeltaR1stJet2ndJet,"Imbalance against deltaR1stJet2ndJet",40,0.4,1.5,50,0,10,"Imbalance","DeltaR");
  hImbalanceEta2ndJet = createTH2(hImbalanceEta2ndJet,"Imbalance against Eta 2nd Jet",40,0.4,1.5,50,0,10,"Imbalance","eta 2nd Jet");
  hImbalanceRatioPhotonPt12 = createTH2(hImbalanceRatioPhotonPt12,"Imbalance against Ratio of first two photons",40,0.4,1.5,100,0,1,"Imbalance","Ratio 12Photon");
  hImbalanceRatioPhotongenPhoton = createTH2(hImbalanceRatioPhotongenPhoton,"Imbalance against Ratio of ph/genPh",40,0.4,1.5,100,0,2,"Imbalance","Ratio Photon/genPhoton");  
  hImbalanceDeltaPhiPhoton2ndJet = createTH2(hImbalanceDeltaPhiPhoton2ndJet,"Imbalance against deltaphit pho 2nd Jet",40,0.4,1.5,100,0,3.2,"Imbalance","delta phi");  
  hImbalancePhoton2Pt = createTH2(hImbalancePhoton2Pt,"Imbalance against photon 2Pt",40,0.4,1.5,100,0,600,"Imbalance","Photon2Pt");  
  hImbalanceDeltaEtaPhoton2ndJet = createTH2(hImbalanceDeltaEtaPhoton2ndJet,"Imbalance against deltaEta pho 2nd Jet",40,0.4,1.5,100,0,5,"Imbalance","delta Eta");  
  hImbalanceNobjPhoton = createTH2(hImbalanceNobjPhoton,"Imbalance against nObjPhoton",40,0.4,1.5,20,0,20,"Imbalance","nObjPhotons");  
  hImbalance2ndJetCorr = createTH2(hImbalance2ndJetCorr,"Imbalance against L1L2L3corr",40,0.4,1.5,100,0,3,"Imbalance","L1L2L3corr"); 
  hImbalance1stGenJetID = createTH2(hImbalance1stGenJetID,"Imbalance against Id of first genJet",40,0.4,1.5,25,0,25,"Imbalance","genJetid"); 
  hImbalance2ndGenJetID = createTH2(hImbalance2ndGenJetID,"Imbalance against Id of second genJet",40,0.4,1.5,25,0,25,"Imbalance","genJetid"); 
  hPhotonEta2ndJet = createTH2(hPhotonEta2ndJet,"Photon Eta against 2nd Jet Eta",40,-1.3,1.3,40,-6,6,"PhotonEta","2ndJetEta"); 
  hImbalanceRatioPhotonPtJetColPt = createTH2(hImbalanceRatioPhotonPtJetColPt,"Imbalance against Ratio of photon pt/JetCol pt",40,0.5,1.5,40,0.6,1.4,"Imbalance","photonPt/ColPt"); 
  
  

  hJet1PtPhoton1Pt2 = createTH2(hJet1PtPhoton1Pt2,"Jet1 Pt and Photon1 Pt 2",100,0,upp_pt,100,130,upp_pt,"p_{T}^{Jet}","p_{T}^{#gamma}");
  hImbalanceJetN2 = createTH2(hImbalanceJetN2,"Imbalance against Jet Multiplicity 2",40,0.4,1.5,100,0,100,"Imbalance","# Jets");
  hImbalancePhoton1Pt2 = createTH2(hImbalancePhoton1Pt2,"Imbalance against Photon1 Pt 2",40,0.4,1.5,50,0,4000,"Imbalance","Photon1 Pt");
  hImbalanceVtx2 = createTH2(hImbalanceVtx2,"Imbalance against number of vertices 2",40,0.4,1.5,50,0,50,"Imbalance","# Vtx");
  hImbalanceJetPt2= createTH2(hImbalanceJetPt2,"Imbalance against jetPt 2",100,0.4,1.5,100,1,3000,"Imbalance","jetPt");
  hImbalanceJetEta2= createTH2(hImbalanceJetEta2,"Imbalance against jet eta 2",100,0.4,1.5,100,0,1.0,"Imbalance","jet Eta");
  hImbalancePhotonEta2 = createTH2(hImbalancePhotonEta2,"Imbalance against #eta_{#gamma} 2",20,0.4,1.5,100,-1.3,1.3,"Imbalance","#eta_{#gamma}");
  hImbalanceWeights2 = createTH2(hImbalanceWeights2,"Imbalance against Weights 2",100,0.4,1.5,200,0,TMath::Power(10,-5),"Imbalance","Weights");
  hImbalanceJet2Pt2 = createTH2(hImbalanceJet2Pt2,"Imbalance against Jet2Pt 2",20,0.4,1.5,500,0,500,"Imbalance","Jet2Pt");
  hImbalanceGenJet2Pt2 = createTH2(hImbalanceGenJet2Pt2,"Imbalance against genJet2Pt 2",20,0.4,1.5,500,0,500,"Imbalance","genJet2Pt");
  hJet2PtGenJet2Pt2 = createTH2(hJet2PtGenJet2Pt2,"Jet2Pt against genJet2Pt 2",500,0,500,500,0,500,"Jet2Pt","genJet2Pt");
  hImbalanceAlpha2 = createTH2(hImbalanceAlpha2,"Imbalance against Alpha 2",20,0.4,1.5,50,0,10,"Imbalance","Alpha");
  hImbalanceDeltaRPhoton1stJet2 = createTH2(hImbalanceDeltaRPhoton1stJet2,"Imbalance against deltaRPhoton1stJet 2",40,0.4,1.5,50,0,10,"Imbalance","DeltaR");
  hImbalanceDeltaRPhoton2ndJet2 = createTH2(hImbalanceDeltaRPhoton2ndJet2,"Imbalance against deltaRPhoton2ndJet 2",40,0.4,1.5,50,0,10,"Imbalance","DeltaR");
  hImbalanceDeltaR1stJet2ndJet2 = createTH2(hImbalanceDeltaR1stJet2ndJet2,"Imbalance against deltaR1stJet2ndJet 2",40,0.4,1.5,50,0,10,"Imbalance","DeltaR");
  hImbalanceEta2ndJet2 = createTH2(hImbalanceEta2ndJet2,"Imbalance against Eta 2nd Jet 2",40,0.4,1.5,50,0,10,"Imbalance","eta 2nd Jet");
  hImbalanceRatioPhotonPt122 = createTH2(hImbalanceRatioPhotonPt122,"Imbalance against Ratio of first two photons2",40,0.4,1.5,100,0,1,"Imbalance","Ratio 12Photon");
  hImbalanceRatioPhotongenPhoton2 = createTH2(hImbalanceRatioPhotongenPhoton2,"Imbalance against Ratio of ph/genPh 2",40,0.4,1.5,100,0,2,"Imbalance","Ratio Photon/genPhoton");  
  hImbalanceDeltaPhiPhoton2ndJet2 = createTH2(hImbalanceDeltaPhiPhoton2ndJet2,"Imbalance against deltaphit pho 2nd Jet 2",40,0.4,1.5,100,0,3.2,"Imbalance","delta phi");
  hImbalancePhoton2Pt2 = createTH2(hImbalancePhoton2Pt2,"Imbalance against photon 2Pt 2",40,0.4,1.5,100,0,600,"Imbalance","Photon2Pt");  
  hImbalanceDeltaEtaPhoton2ndJet2 = createTH2(hImbalanceDeltaEtaPhoton2ndJet2,"Imbalance against deltaEta pho 2nd Jet 2",40,0.4,1.5,100,0,5,"Imbalance","delta Eta");  
  hImbalanceNobjPhoton2 = createTH2(hImbalanceNobjPhoton2,"Imbalance against nObjPhoton 2",40,0.4,1.5,20,0,20,"Imbalance","nObjPhotons");  
  hImbalance2ndJetCorr2 = createTH2(hImbalance2ndJetCorr2,"Imbalance against L1L2L3corr 2",40,0.4,1.5,100,0,3,"Imbalance","L1L2L3corr");  
  hImbalance1stGenJetID2 = createTH2(hImbalance1stGenJetID2,"Imbalance against Id of first genJet 2",40,0.4,1.5,25,0,25,"Imbalance","genJetid"); 
  hImbalance2ndGenJetID2 = createTH2(hImbalance2ndGenJetID2,"Imbalance against Id of second genJet 2",40,0.4,1.5,25,0,25,"Imbalance","genJetid"); 
  hPhotonEta2ndJet2 = createTH2(hPhotonEta2ndJet2,"Photon Eta against 2nd Jet Eta2",40,-1.3,1.3,40,-6,6,"PhotonEta","2ndJetEta"); 
  hImbalanceRatioPhotonPtJetColPt2 = createTH2(hImbalanceRatioPhotonPtJetColPt2,"Imbalance against Ratio of photon pt/JetCol pt 2",40,0.5,1.5,40,0.6,1.4,"Imbalance","photonPt/ColPt"); 

  // 2d: Histogram for jet PT in small response region
  hSmallJetPtResponse = createTH2(hSmallJetPtResponse,"Jet Pt in low response region against Response",100,0,upp_pt,100,0,1,"p_{T,1} (GeV)","Response");

  // 2d: Histogram for jet PT in small response region
  hPhotonEta_high_pt_reg_Response = createTH2(hPhotonEta_high_pt_reg_Response,"Photon Eta in low response region against Response",48,-3,3,10,0,1,"#eta^{#gamma}","Response");

  // 2d: Histogram for Response against EMF for leading Jet
  hResponseJetEMF = createTH2(hResponseJetEMF,"Response against EMF (leading Jet)",50,0,1,50,0,1,"Response","EMF");
 
  // 2d: Histogram for Response against FHPD for leading Jet
  hResponseJetFHPD = createTH2(hResponseJetFHPD,"Response against FHPD (for leading Jet)",50,0,1,50,0,1,"Response","FHPD");

  // 2d: Histogram for Response against FRBX for leading Jet
  hResponseJetFRBX = createTH2(hResponseJetFRBX,"Response against FRBX (leading Jet)",50,0,1,50,0,1,"Response","FRBX");

  // 2d: Histogram for Response against EMF for Photon
  hResponsePhotonEMF = createTH2(hResponsePhotonEMF,"Response against EMF (Photon)",50,0,1,50,0,1,"Response","EMF");

  hResponseLeadJet = createTH2(hResponseLeadJet,"Response against lead_jet Number",50,0,2,100,0,100,"Response","lead_jet");

  // 2d: Histogram for Response against FHPD for Photon
  hResponsePhotonFHPD = createTH2(hResponsePhotonFHPD,"Response against FHPD (Photon)",50,0,1,50,0,1,"Response","FHPD");

  // 2d: Histogram for Response against FRBX for Photon
  hResponsePhotonFRBX = createTH2(hResponsePhotonFRBX,"Response against FRBX (Photon)",50,0,1,50,0,1,"Response","FRBX");

   // 2d: Histogram for Response against FRBX for Photon
  hResponsePhotonEta = createTH2(hResponsePhotonEta,"Response against Photon Eta)",40,0.4,1.5,60,-3,-3,"Response","Photon Eta");

  // 2d: Histogram for Response against Pt of Photon from Photon collection/pt from photon from Jet collection
  hResponsePhotonPtRatio = createTH2(hResponsePhotonPtRatio,"Response against Photon pt ratio jet and photon collection",50,0,1,50,0,2,"Response","Pt ratio");
   
  // 2d: Pt of first two reconstructed Photons
  hPhoton12Pt = createTH2(hPhoton12Pt,"P_{t} of first two photons",100,0,upp_pt,100,0,upp_pt,"p_{T,1} GeV","p_{T,2} GeV");

  // 2d: Eta against Pt of first Jet
  hJetEtaPt1 = createTH2(hJetEtaPt1,"Eta and P_{t} of first Jet",100,-5.5,5,100,0,upp_pt,"Eta Jet","p_{T}^{Jet_{1}} GeV");

  

 
 
  
  
 

  

  

  

  // 2d: Response against Pt of first Photon before deltaEta cut
  hResponsePt = createTH2(hResponsePt,"Response and P_{t} of first Photon before #Delta #eta cut ",100,0,2.5,100,0,3000,"Response","p^{#gamma}_{T,1} GeV");

  // 2d: Response against Pt of first Photon for deltaEta < 0.5
  hResponsePt_0_5_Deta = createTH2(hResponsePt_0_5_Deta,"Response and P_{t} of first Photon for #Delta #eta < 0.5",100,0,2.5,100,0,upp_pt,"Response","p^{#gamma}_{T,1} GeV");
  
  // 2d: Response against Pt of first Photon for deltaEta < 1.1
  hResponsePt_1_1_Deta = createTH2(hResponsePt_1_1_Deta,"Response and P_{t} of first Photon for #Delta #eta < 1.1",100,0,2.5,100,0,upp_pt,"Response","p^{#gamma}_{T,1} GeV");

  // 2d: Response against Pt of first Photon for deltaEta < 1.5
  hResponsePt_1_5_Deta = createTH2(hResponsePt_1_5_Deta,"Response and P_{t} of first Photon for #Delta #eta < 1.5",100,0,2.5,100,0,upp_pt,"Response","p^{#gamma}_{T,1} GeV");
 
  // 2d: Response against Pt of first Photon before ecal cut
  hResponsePtBeforeEcal = createTH2(hResponsePtBeforeEcal,"Response and P_{t} of first Photon before Ecal cut",100,0,2.5,100,0,upp_pt,"Response","p^{#gamma}_{T,1} GeV");

  // 2d: Response against Jet Energy Fraction (scalar sum)
  hResponseEnergyFraction = createTH2(hResponseEnergyFraction,"Response against #Sigma p^{Jet}_{2,3,4,5}/p^{Jet}_{1} (scalar sum)",100,0,2.5,50,0,5,"Response","#frac{#Sigma p^{Jet}_{2,3,4,5}}{p^{Jet}_{1}}");
  
  // 2d: Response against Jet Energy Fraction (vectorial sum)
  hResponseVecEnergyFraction = createTH2(hResponseVecEnergyFraction,"Response against (#Sigma p^{Jet}_{2,3,4,5})/p^{Jet}_{1} (vectorial sum)",100,0,2.5,50,0,5,"Response","#frac{#Sigma p^{Jet}_{2,3,4,5}}{p^{Jet}_{1}}");

  // 2d: Response against L1-correction
  hResponseL1 = createTH2(hResponseL1,"Response against L1 correction",100,0,2.5,32,0.4,1.2,"Response","L1 correction");
  
  // 2d: Response against L2L3-correction
  hResponseL2L3= createTH2(hResponseL2L3,"Response against L2L3 - correction",100,0,2.5,32,1,1.4,"Response","L2L3 - correction");

  // Photon Pt against Vtx
  hPhotonPtVtx= createTH2(hPhotonPtVtx,"Photon Pt against Number of Vertices",100,0,upp_pt,60,0,60,"Photon PT","NVtx");

  // Photon Pt against Vtx
  
  for(int i=0; i<numTrigger-1;i++){
    if(i<numTrigger-2) sprintf(tname,"Photon Pt against Number of Vertices for %4.1f GeV < p_{T} < %4.1f GeV ",ptBins[i],ptBins[i+1]);
    else    sprintf(tname,"Photon Pt against Number of Vertices for %4.1f GeV < p_{T}",ptBins[i]);
    hPhotonPtVtxBinned[i]= createTH2(hPhotonPtVtxBinned[i],tname,100,0,upp_pt,60,0,60,"Photon PT","NVtx");
  }
  
  //PU distribution versus Photon Pt before all cuts

  hPhotonPtVtxbeforeCUTS=createTH2(hPhotonPtVtx,"Photon Pt against Number of Vertices (before CUTS)",1200,1,upp_pt,60,0,60,"Photon PT","NVtx");

}
