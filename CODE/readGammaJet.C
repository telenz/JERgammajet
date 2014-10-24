#include "myDeclarations.h"
#include "myClasses.h"

#include <iostream>
#include "TChain.h"
#include "TString.h"


int readGammaJet(int nEvents) {

  TString DataFilename, PUVersion;
  chain  = new TChain("GammaJetTree");

  // Determine from which root file to read;
  if(date == 2012){
    if(isMC || testClosure || QCDUncertaintyEvaluation){
      if(set == 1){
	DataFilename = "/scratch/hh/lustre/cms/user/telenz/mc/PhotonJetTuple2012/allTogether/OnlyTightPhotons/";
	if(jetType==1)      DataFilename += "ak5FastPF_*.root";  
	else if(jetType==2) DataFilename += "ak5PFCHS_*.root";  
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
      else if(set == 2 || set == 3 || set == 4){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/pythia_flat_53/MCPhoton2012/OnlyTightPhotons/";
	if(jetType==1)      DataFilename += "ak5FastPF.root";
	else if(jetType==2) DataFilename += "ak5PFCHS.root";    
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
      else if(set == 5){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/pythia_flat_535_START53_V20_WithHasPixelSeed/MCPhoton2012/OnlyTightPhotons/";
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/pythia_flat_535_START53_V22/MCPhoton2012/OnlyTightPhotons/";
	DataFilename = "/nfs/dust/cms/user/tlenz/mc/pythia_flat_535_START53_V22/MCPhoton2012/OnlyTightPhotons/";
	if(jetType==1)      DataFilename += "ak5FastPF*.root";
	else if(jetType==2) DataFilename += "ak5PFCHS*.root"; 
	else if(jetType==4) DataFilename += "ak7PFCHS*.root";    
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }

    }
    else if(!isMC){
      if(set == 1){
	DataFilename = "/scratch/hh/lustre/cms/user/telenz/data/PhotonJetTuple2012AB/";      
	if(jetType==1)      DataFilename += "ak5FastPF_*.root";
	else if(jetType==2) DataFilename += "ak5PFCHS_*.root"; 
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
      else if(set == 2){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPF*.root";
	else if(jetType==2) DataFilename += "ak5PFCHS*.root";     
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
      else if(set == 3){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC/otherSamples/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPFAB.root";
	else if(jetType==2) DataFilename += "ak5PFCHSAB.root";     
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
      else if(set == 4){
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC/OnlyTightPhotons/";      
	if(jetType==1)      DataFilename += "ak5FastPFC*.root";
	else if(jetType==2) DataFilename += "ak5PFCHSC*.root";     
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
      else if(set == 5){
	//DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2012ABC_5_3_5_FT_53_V21/OnlyTightPhotons/";     
	DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/2012/PhotonJetTuple2012ABCDrereco_5_3_5_FT_53_V21_AN4/OnlyTightPhotons/";
	DataFilename = "/nfs/dust/cms/user/tlenz/data/2012/PhotonJetTuple2012ABCDrereco_5_3_5_FT_53_V21_AN4/OnlyTightPhotons/";
	if(jetType==1)      DataFilename += "ak5FastPF_*.root";
	else if(jetType==2) DataFilename += "ak5PFCHS_*.root";     
	else if(jetType==4) DataFilename += "ak7PFCHS_*.root";     
	else{
	  cout<<"No input File available!"<<endl;
	  return 0;
	}
      }
    }
    else cout<<"No such \"type\" available. Please choose a different number for variable type."<<endl;
  }
  else if (date == 2011){
    if(!isMC){
      DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/data/PhotonJetTuple2011/OnlyTightPhotons/";
      if(jetType==1) DataFilename += "ak5FastJet.root";
      else{
	cout<<"No input File available!"<<endl;
	return 0;
      }
    }
    else if(isMC){
      DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/mc_2011/OnlyTightPhotons/";
      if(jetType==1) DataFilename += "ak5FastJet.root";
      else{
	cout<<"No input File available!"<<endl;
	return 0;
      }
    }
  }
  else  cout<<"No input File available! (no none 2012 data available)"<<endl;
  cout<<"filename = "<<DataFilename<<endl; 
  chain->Add(DataFilename); 

  if((QCDUncertaintyEvaluation && !isMC) || addQCD){
    //if(jetType==2) chain->Add("/scratch/hh/dust/naf/cms/user/telenz/mc/QCDenriched_pythia_PtBinned_START53_V20/OnlyTightPhotons/ak5PFCHS*.root"); 
    if(jetType==2) chain->Add("/scratch/hh/dust/naf/cms/user/telenz/mc/QCDenriched_pythia_PtBinned_START53_V22/OnlyTightPhotons/ak5PFCHS*.root"); 
    else{
      cout<<"No input File available!"<<endl;
      return 0;
    }
    
    cout<<endl<<"QCD Sample was Added!!!!!!!"<<endl<<endl;
  }

  if(addWJet && isMC){
    //if(jetType==1) chain->Add("/scratch/hh/dust/naf/cms/user/telenz/mc/qcd_pythia_flat_535_START53_V20/OnlyTightPhotons/ak5FastPF*.root"); 
    if(jetType==2) chain->Add("/scratch/hh/dust/naf/cms/user/telenz/mc/WJets_madgraph_HTbinned_START53_V20/OnlyTightPhotons/ak5PFCHS*.root"); 
    if(jetType==4) chain->Add("/scratch/hh/dust/naf/cms/user/telenz/mc/WJets_madgraph_HTbinned_START53_V20/OnlyTightPhotons/ak7PFCHS*.root"); 
    cout<<endl<<"W+Jet Sample was Added!!!!!!!"<<endl<<endl;
  }

  // Set branch addresses
  chain->SetBranchAddress("CrossSection",&crossSection);
  chain->SetBranchAddress("NobjPhoton",&nobjPhoton);
  chain->SetBranchAddress("NobjGenPhoton",&nobjGenPhoton);
  chain->SetBranchAddress("PhotonPt",photonPt);
  chain->SetBranchAddress("PhotonE",photonE);
  chain->SetBranchAddress("PhotonEta",photonEta);
  chain->SetBranchAddress("PhotonPhi",photonPhi);
  chain->SetBranchAddress("GenPhotonPt",genPhotonPt);
  chain->SetBranchAddress("GenPhotonEta",genPhotonEta);
  chain->SetBranchAddress("GenPhotonPhi",genPhotonPhi);
  chain->SetBranchAddress("PhotonIDTight",tight);
  chain->SetBranchAddress("PhotonIDLoose",loose);

  chain->SetBranchAddress("NobjJet",&nobjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetE",jetE);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  chain->SetBranchAddress("JetCorrL1",jetCorrL1);
  if(date == 2012) chain->SetBranchAddress("JetCorrUncert",jetUncert);
  chain->SetBranchAddress("JetPhi",jetPhi);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetIDTight",jetIDTight);
  chain->SetBranchAddress("JetIDLoose",jetIDLoose);
  chain->SetBranchAddress("VtxN",&vtxN);

  chain->SetBranchAddress("PassesECALDeadCellBEFilter",&EcalBEFilter);
  chain->SetBranchAddress("PassesECALDeadCellTPFilter",&EcalTPFilter);

  chain->SetBranchAddress("RunNumber",&runNum);
  chain->SetBranchAddress("LumiBlockNumber",&lumiNum);
  chain->SetBranchAddress("EventNumber",&eventNum);  

  chain->SetBranchAddress("PhotonIsoECAL04",photonIsoEcal);  
  chain->SetBranchAddress("PhotonIsoHCAL04",photonIsoHcal);
  chain->SetBranchAddress("PhotonIsoTrk04",photonIsoTrk);
  chain->SetBranchAddress("PhotonSigmaIetaIeta",photonSigmaIetaIeta);
  chain->SetBranchAddress("PhotonHadronicOverEM",photonHadronicOverEM);
  chain->SetBranchAddress("PhotonHasPixelSeed",photonHasPixelSeed);

  chain->SetBranchAddress("JetFRBX",jetFRBX);
  chain->SetBranchAddress("JetFHPD",jetFHPD);
  chain->SetBranchAddress("JetEMF",jetEMF);
  
  if(isMC || testClosure || QCDUncertaintyEvaluation)      chain->SetBranchAddress("NewWeight",&weight);
  else if(!isMC || date == 2011) chain->SetBranchAddress("Weight",&weight);

  chain->SetBranchAddress("Rho",&rho);
  chain->SetBranchAddress("PUMCNumTruth",&PUMCNumTruth);

  if(isMC){
    chain->SetBranchAddress("Met_T1",&met_T1);
    chain->SetBranchAddress("MetPhi_T1",&metPhi_T1);
  }
  else if(!isMC){
    chain->SetBranchAddress("Met_T1R",&met_T1);
    chain->SetBranchAddress("MetPhi_T1R",&metPhi_T1);
  }

  chain->SetBranchAddress("GenJetColPt",genJetPt);
  chain->SetBranchAddress("GenJetColEta",genJetEta);
  chain->SetBranchAddress("GenJetColPhi",genJetPhi);
  chain->SetBranchAddress("GenPartId_algo",genJetID_algo);
  chain->SetBranchAddress("GenPartId_phys",genJetID_phys);
  chain->SetBranchAddress("GenJetColJetIdx",genJetColJetIdx);
  chain->SetBranchAddress("NobjGenJet",&nobjGenJet);
  chain->SetBranchAddress("GenEvtScale",&genEvtScale);

  chain->SetBranchAddress("HltPhoton20iso",&hltPhoton[0]);
  chain->SetBranchAddress("HltPhoton30iso",&hltPhoton[1]);
  chain->SetBranchAddress("HltPhoton50iso",&hltPhoton[2]);
  chain->SetBranchAddress("HltPhoton75iso",&hltPhoton[3]);
  chain->SetBranchAddress("HltPhoton90iso",&hltPhoton[4]);  

  if(date == 2012){
    chain->SetBranchAddress("HltPhoton135",&hltPhoton[5]);
    chain->SetBranchAddress("HltPhoton150",&hltPhoton[6]);
    chain->SetBranchAddress("HltPhoton160",&hltPhoton[7]);
  }

  if(isMC){
    if(PUreweighting){

      if(set == 1 && date == 2012)      PUVersion = "PUWeights/PUData_2012/weightsCMSSW5_2_5/";
      else if(set == 2 && date == 2012) PUVersion = "PUWeights/PUData_2012/weightsCMSSW5_3_3/";
      else if(set == 3 && date == 2012) PUVersion = "PUWeights/PUData_2012/weightsCMSSW5_3_3_AB/";
      else if(set == 4 && date == 2012) PUVersion = "PUWeights/PUData_2012/weightsCMSSW5_3_3_C/";
      //else if(set == 5 && date == 2012 && !addQCD) PUVersion = "PUData_2012ABCD/weightsCMSSW5_3_5/";
      //else if(set == 5 && date == 2012 && !addQCD) PUVersion = "PUData_2012ABCrereco/weightsCMSSW_5_3_5/";
      //else if(set == 5 && date == 2012 && !addQCD) PUVersion = "PUData_2012ABCrereco/weightsCMSSW_5_3_5/";
      //else if(set == 5 && date == 2012 && !addQCD) PUVersion = "PUData_2012ABCD/weightsCMSSW5_3_5_MBXS69300_TEST/";
      else if(set == 5 && date == 2012 && !addQCD) PUVersion = "PUWeights/PUData_2012ABCDrereco/cmssw5_3_5_Tag04_02_02LumiCalc2FineBinning_MBXS69400/weights_60Bins/";
      else if(set == 5 && date == 2012 && addQCD)  PUVersion = "PUWeights/PUData_2012ABCD/weightsCMSSW_5_3_5_plusQCD/";
      //else if(set == 5 && date == 2012) PUVersion = "PUData_2012ABCD/weightsCMSSW5_3_5_WithWeights/";
      //else if(set == 5 && date == 2012) PUVersion = "PUData_2012ABCD/weightsNVTX/";
      else if(date == 2011) PUVersion = "PUData_2011/newVersion/";
      else cout<<"No PU data histograms available !!! (look in CODE/readGammaJets.C)"<<endl;

      cout<<endl<<"PU Weighst taken from following folder:  "<<PUVersion<<endl<<endl;
 
  
      // Read PU Weight-distribution  (seperately for different Triggers)
      hPUWeight[0] =  readTH1(PUVersion + "PUWeights_0.root","PUWeight0","PUWeight0");
      hPUWeight[1] =  readTH1(PUVersion + "PUWeights_1.root","PUWeight1","PUWeight1");
      hPUWeight[2] =  readTH1(PUVersion + "PUWeights_2.root","PUWeight2","PUWeight2");
      hPUWeight[3] =  readTH1(PUVersion + "PUWeights_3.root","PUWeight3","PUWeight3");
      hPUWeight[4] =  readTH1(PUVersion + "PUWeights_4.root","PUWeight4","PUWeight4");
  
      if(date==2012){
	hPUWeight[5] =  readTH1(PUVersion + "PUWeights_5.root","PUWeight5","PUWeight5");
	hPUWeight[6] =  readTH1(PUVersion + "PUWeights_6.root","PUWeight6","PUWeight6");
	hPUWeight[7] =  readTH1(PUVersion + "PUWeights_7.root","PUWeight7","PUWeight7");
      }
    }
    else cout<<endl<<"No PU reweighting applied!!!"<<endl<<endl;
  }


  // If tree contains less entries than nEvents use less events
  nMax = nEvents;
  if( chain->GetEntries() < nEvents ) nMax = chain->GetEntries();
  cout<<"There are "<<chain->GetEntries()<<" events in this sample!"<<endl;

  return 1;
  
}
