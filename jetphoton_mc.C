//  Script to calculate JET Energy Response in Gamma Jet Events --
//  (see: CMS AN-2010/076 and CMS AN-2010/421)                  --
// ---------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"


const bool isMC    = true;

//class CsysUncerMC;

#include "CODE/myDeclarations.h"
#include "CODE/myClasses.h"
#include "CODE/myFunctions.h"
#include "CODE/myClasses.C"
#include "CODE/myFunctions.C"
#include "CODE/calcPath.C"
#include "CODE/bookHistos.C"
#include "CODE/readGammaJet.C"

#include "CODE/controlPlots.C"
#include "CODE/applyCuts.C"
#include "CODE/extrapolate.C"



void calcSample();
void calcScale();

  
// ---- Function implementations ---------------
// ---------------------------------------------

// Run main script
// This is the actual command to be typed
// used in the ROOT session.
// ---------------------------------------------
int jetphoton_mc(int nEvts, int step) {

  
  // -------------------------------------
  if(date == 2011 && numTrigger!=5){
    cout<<"Please correct the number of Triggers used!"<<endl;  
    return 0;
  }
  if(date == 2012 && numTrigger!=7){
    cout<<"Please correct the number of Triggers used!"<<endl;  
    return 0;
  }

  //-------- print out of Information --------------
  cout<<endl;
  if(!isMC)     cout<<"Data!!"<<endl<<endl;
  else if(isMC) cout<<"MC!!"<<endl<<endl;

  if(date == 2012) cout<<"date = 2012"<<endl<<endl;
  else if(date == 2011) cout<<"date = 2011"<<endl<<endl;
  
  if(jetType == 1)      cout<<"jetType = ak5PF Jets "<<endl<<endl<<endl;
  else if(jetType == 2) cout<<"jetType = ak5PFCHS Jets"<<endl<<endl<<endl;
  else if(jetType == 4) cout<<"jetType = ak7PFCHS Jets"<<endl<<endl<<endl;

  
  calcPath(step); 
  bookHistos();
  readGammaJet(nEvts);
  time_t start,end;
  time (&start);

  if(step==1){
    cout<<"calcSample() is executed!"<<endl<<endl;
    calcSample();
  }
  if(step==2){
    cout<<"calcScale() is executed!"<<endl<<endl;

    if(detJER == 1)      cout<<endl<<"Method is : Gaus 2 sigma!!"<<endl;
    else if(detJER == 2) cout<<endl<<"Method is : RMS 95!!"<<endl;
    else if(detJER == 3) cout<<endl<<"Method is : RMS 99!!"<<endl;

    calcScale();
  }
  //----------------
  time (&end);
  
  cout<<"elapsed time:"<<difftime(end,start)/60<<" minutes"<<endl;

  return 0;
}


// Aplly all cuts and fill histograms
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void calcSample() {

  
  initializeControlPlots();
 
  for(int j=0; j<nEtaBins; j++){
    for(int i=0; i<nPtBins; i++){
      for(int k=0; k<nAlphaBins; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 	
	JetImbalanceJetHemisphere[i][j][k]        = new CResponse(3);
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 
	JetImbalancePhotonHemisphere[i][j][k]     = new CResponse(3);
	JetIntrinsic[i][j][k]                     = new CResponse(2);
	JetIntrinsicle5[i][j][k]                  = new CResponse(2);
	JetIntrinsicle10gt5[i][j][k]              = new CResponse(2);
	JetIntrinsicle15gt10[i][j][k]             = new CResponse(2);
	JetIntrinsicle20gt15[i][j][k]             = new CResponse(2);
	JetIntrinsicgt20[i][j][k]                 = new CResponse(2);	
   
      } 
    }  
  }

  //------------------------------------------------------------------------------------------------------------
  // Loop over nMax entries, read variables and fill histogram

  std::cout << "Processing events" << std::endl;
  
  
  for(int n = 0; n < nMax; n++) {
    
    if( (n+1)%1000000 == 0 ) std::cout << "Event " << (n+1) << std::endl;
   
    // Get this event i.e. entry n and fill values into variables
    chain->GetEntry(n);
      
    // Calculating PUWeight for this event   
    if(PUreweighting && !testClosure){
      for(int i=0; i<numTrigger-1; i++){
	if(photonPt[0] >= ptBins[i] && photonPt[0] < ptBins[i+1]){
	
	  PUWeight = hPUWeight[i]->GetBinContent(hPUWeight[i]->FindBin(PUMCNumTruth));
	  hPUgenMC[i]     -> Fill(PUMCNumTruth,weight*PUWeight);
	  break;
	}
	else if(photonPt[0]>=ptBins[numTrigger-1]){
	
	  PUWeight = hPUWeight[numTrigger-1]->GetBinContent(hPUWeight[numTrigger-1]->FindBin(PUMCNumTruth));
	  hPUgenMC[numTrigger-1]    -> Fill(PUMCNumTruth,weight*PUWeight);
	  break;
	}
      }
    }
    else PUWeight=1;
      
    
    //---------------------------------------------------------------------------------------------
    // Sort Jets and apply all cuts with following function
    bool testVar = applyCuts();
    if(!testVar) continue;    
    //---------------------------------------------------------------------------------------------
       
    // Fill Response functions for whole sample    
    for(int k=0; k<nAlphaBins; k++){

      if(alpha >= alphaBins[k] && alpha < alphaBins[k+1]){ 
		
	for(int j=0; j<nEtaBins; j++){
	  
	  if(std::abs(jetEta[corrJets.idx(idx1stJet)])>= etaBins[j] && std::abs(jetEta[corrJets.idx(idx1stJet)]) < etaBins[j+1] ){
	    
	    
	    for(int i=0; i<nPtBins; i++){ 
	      
	      if(photonPt[0] >= ptBins[i] && photonPt[0] < ptBins[i+1]){
		//if(sumVector.Pt() >= ptBins[i] && sumVector.Pt() < ptBins[i+1]){
		//if(genJetPt[gen1stJetidx] >= ptBins[i] && genJetPt[gen1stJetidx] < ptBins[i+1]){
		//if(jetPt[corrJets.idx(idx1stJet)] >= ptBins[i] && jetPt[corrJets.idx(idx1stJet)] < ptBins[i+1]){
	   	
		float deltaphi2ndJet1stJet = 0;
		float deltaphi2ndJetPhoton = 0;
		if(idx2ndJet != -1){
		 
		  deltaphi2ndJet1stJet = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(idx2ndJet)]-jetPhi[corrJets.idx(idx1stJet)]));
		  deltaphi2ndJetPhoton = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(idx2ndJet)]-photonPhi[0])); 
		
	       	}
		if(gen1stJetidx < 0) count1stJetNotMatched += 1;
					 		
		if(deltaphi2ndJetPhoton > deltaphi2ndJet1stJet){

		  JetResponseJetHemisphere[i][j][k] -> hResponse      -> Fill(response,weight*PUWeight);
		  JetResponseJetHemisphere[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);
		  JetResponseJetHemisphere[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);	
		  JetResponseJetHemisphere[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(idx1stJet)],weight*PUWeight);
		  JetResponseJetHemisphere[i][j][k] -> hPtLeadingPhoton -> Fill(photonPt[0],weight*PUWeight);
		 

		  if(gen1stJetidx >= 0){
		    JetResponseJetHemisphere[i][j][k]  -> hgenPtLeadingJet  -> Fill(genJetPt[gen1stJetidx],weight*PUWeight);  
		    JetImbalanceJetHemisphere[i][j][k] -> hResponse         -> Fill(imbalance,weight*PUWeight);
		    JetImbalanceJetHemisphere[i][j][0] -> hPt               -> Fill(photonPt[0],weight*PUWeight);
		    JetImbalanceJetHemisphere[i][j][k] -> hAlpha            -> Fill(alpha,weight*PUWeight);
		  }		  		  	
		}  
		else{
		  JetResponsePhotonHemisphere[i][j][k] -> hResponse      -> Fill(response,weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hPtLeadingJet  -> Fill(jetPt[corrJets.idx(idx1stJet)],weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hPtLeadingPhoton -> Fill(photonPt[0],weight*PUWeight);

		  if(gen1stJetidx >= 0){
		    JetResponsePhotonHemisphere[i][j][k] -> hgenPtLeadingJet  -> Fill(genJetPt[gen1stJetidx],weight*PUWeight);  
		    JetImbalancePhotonHemisphere[i][j][k]-> hResponse         -> Fill(imbalance,weight*PUWeight); 
		    JetImbalancePhotonHemisphere[i][j][0]-> hPt               -> Fill(photonPt[0],weight*PUWeight); 
		    JetImbalancePhotonHemisphere[i][j][k]-> hAlpha            -> Fill(alpha,weight*PUWeight);
		  } 
		}

		if(gen1stJetidx >= 0){
		  JetIntrinsic[i][j][k] -> hResponse         -> Fill(intrinsic,weight*PUWeight); 
		  JetIntrinsic[i][j][0] -> hPt               -> Fill(photonPt[0],weight*PUWeight);	
		  JetIntrinsic[i][j][k] -> hAlpha            -> Fill(alpha,weight*PUWeight);
		  JetIntrinsic[i][j][k] -> hPtLeadingJet     -> Fill(jetPt[corrJets.idx(idx1stJet)],weight*PUWeight);
		  JetIntrinsic[i][j][k] -> hgenPtLeadingJet  -> Fill(genJetPt[gen1stJetidx],weight*PUWeight);  
		  JetIntrinsic[i][j][k] -> hPtLeadingPhoton  -> Fill(photonPt[0],weight*PUWeight); 
		}
	
		hNPU -> Fill(PUMCNumTruth, weight*PUWeight);
		hNPV -> Fill(vtxN, weight*PUWeight);	
		nocut = nocut +1; 
		break;
	      }
	    }
	    break;
	  }
	}
	break;
      }
    }
       
  } // End of loop over entries

  TString TotFilename;  

  // Scale Response functions to MC statistics
  for(int i=0;i<nPtBins;i++){
    for(int j=0; j<nEtaBins; j++){
      
      // Photon Pt Histograms
      saveObject(JetResponseJetHemisphere[i][j][0]->hPt, RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo");   
      saveObject(JetResponsePhotonHemisphere[i][j][0]->hPt, RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo"); 
      saveObject(JetImbalanceJetHemisphere[i][j][0]->hPt , RootPath + "hPt_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo"); 
      saveObject(JetImbalancePhotonHemisphere[i][j][0]->hPt, RootPath + "hPt_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo"); 
      saveObject(JetIntrinsic[i][j][0]->hPt, RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo"); 
      
            
 
      for(int k=0; k<nAlphaBins; k++){	

	// Alpha Histograms
	saveObject(JetResponseJetHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_jet_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root","histo");
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_photon_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root","histo");
	saveObject(JetImbalanceJetHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_imbalance_jet_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root", "histo");
	saveObject(JetImbalancePhotonHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_imbalance_photon_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root", "histo");
	saveObject(JetIntrinsic[i][j][k]->hAlpha,RootPath+"hAlpha_intrinsic_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root", "histo");

	// Response Histograms
	saveObject(JetResponseJetHemisphere[i][j][k]->hResponse, RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	saveObject(JetImbalanceJetHemisphere[i][j][k]->hResponse, RootPath + "response_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hResponse, RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	saveObject(JetImbalancePhotonHemisphere[i][j][k]->hResponse, RootPath + "response_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	saveObject(JetIntrinsic[i][j][k]->hResponse, RootPath + "response_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");

	// Leading Jet Pt
	saveObject(JetResponseJetHemisphere[i][j][k]->hPtLeadingJet, RootPath + "hPtLeadingJet_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");   
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hPtLeadingJet, RootPath + "hPtLeadingJet_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo"); 
	saveObject(JetIntrinsic[i][j][k]->hPtLeadingJet, RootPath + "hPtLeadingJet_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo"); 
	saveObject(JetResponseJetHemisphere[i][j][k]->hgenPtLeadingJet, RootPath + "hgenPtLeadingJet_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");   
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hgenPtLeadingJet, RootPath + "hgenPtLeadingJet_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo"); 
	saveObject(JetIntrinsic[i][j][k]->hgenPtLeadingJet, RootPath + "hgenPtLeadingJet_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo"); 

	// Leading Photon Pt
	saveObject(JetResponseJetHemisphere[i][j][k]->hPtLeadingPhoton, RootPath + "hPtLeadingPhoton_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");   
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hPtLeadingPhoton, RootPath + "hPtLeadingPhoton_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo"); 
	saveObject(JetIntrinsic[i][j][k]->hPtLeadingPhoton, RootPath + "hPtLeadingPhoton_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo"); 
	saveObject(JetIntrinsic[i][j][k]->hPtLeadingPhoton, RootPath + "hPtLeadingPhoton_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");

      }
    }
  }

  // Save some control plots
  plotControlPlots();


  delete chain;
  filestr.close();   
  //EventVariables.close();

  // Special files to write strange events in extra .txt files
  cutflow.open(RootPath + "cutflow.txt");

  cout<<"CUTFLOW (MC):"<<endl<<"cut1 = "<<cut1<<endl<<"cut2 = "<<cut2<<endl<<"cut3 = "<<cut3<<endl<<"cut4 = "<<cut4<<endl<<"cut5 = "<<cut5<<endl<<"cut6 = "<<cut6<<endl<<"cut7 = "<<cut7<<endl<<"cut8 = "<<cut8<<endl<<"cut9 = "<<cut9<<endl<<"cut10 = "<<cut10<<endl<<"cut11 = "<<cut11<<endl;
  cutflow<<"CUTFLOW:"<<endl<<"cut1 = "<<cut1<<endl<<"cut2 = "<<cut2<<endl<<"cut3 = "<<cut3<<endl<<"cut4 = "<<cut4<<endl<<"cut5 = "<<cut5<<endl<<"cut6 = "<<cut6<<endl<<"cut7 = "<<cut7<<endl<<"cut8 = "<<cut8<<endl<<"cut9 = "<<cut9<<endl<<"cut10 = "<<cut10<<endl<<"cut11 = "<<cut11<<endl;
  for(int i=0; i<numTrigger; i++){
    cout<<"cut12["<<i<<"] = "<<cut12[i]<<endl;
    cutflow<<"cut12["<<i<<"] = "<<cut12[i]<<endl;
  }
  for(int i=0; i<3; i++){
    cout<<"cut13["<<i<<"] = "<<cut13[i]<<endl;
    cutflow<<"cut13["<<i<<"] = "<<cut13[i]<<endl;
  }
  cutflow<<"cut14 = "<<cut14<<endl<<"cut15 = "<<cut15<<endl;
  cout<<"cut14 = "<<cut14<<endl<<"cut15 = "<<cut15<<endl;
  cout<<"nocut = "<<nocut<<endl;
  cutflow<<"nocut = "<<nocut<<endl;

  cout<<"number of all events = "<<cut1+cut2+cut3+cut4+cut5+cut6+cut7+cut8+cut9+cut10+cut11+cut12[0]+cut12[1]+cut12[2]+cut12[3]+cut12[4]+cut12[5]+cut12[6]+cut13[0]+cut13[1]+cut13[2]+cut14+cut15+nocut<<endl<<endl;
  cutflow<<"number of all events = "<<cut1+cut2+cut3+cut4+cut5+cut6+cut7+cut8+cut9+cut10+cut11+cut12[0]+cut12[1]+cut12[2]+cut12[3]+cut12[4]+cut12[5]+cut12[6]+cut13[0]+cut13[1]+cut13[2]+cut14+cut15+nocut<<endl<<endl;
 
  cout<<endl<<"count1stJetNotMatched = "<<count1stJetNotMatched<<endl<<endl;
  cout<<"no2ndJetinEvent = "<<no2ndJetinEvent<<endl<<endl;
  cutflow.close();
}

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// Now take the response functions and fit to the core region a gauss function 
// so that one gets the jet energy scale and jet energy resolution


// Calculate sigma and mean for all Response functions
void calcScale(){


  initializeControlPlotsInCalcScale();

  std::cout << "Read Response Histograms!" << std::endl<< std::endl;

  for(int j=0; j<nEtaBins; j++){
    for(int i=0; i<nPtBins; i++){
      for(int k=0; k<nAlphaBins; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 	
	JetImbalanceJetHemisphere[i][j][k]        = new CResponse(3);
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 
	JetImbalancePhotonHemisphere[i][j][k]     = new CResponse(3);
	JetIntrinsic[i][j][k]                     = new CResponse(2);
	
      }
    }
  }

  TString TotFilename;  
  for(int i=0;i<nPtBins;i++){
    for(int j=0; j<nEtaBins; j++){
      
      JetResponseJetHemisphere[i][j][0]->hPt = readTH1(RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo","histo");
      JetResponsePhotonHemisphere[i][j][0]->hPt = readTH1(RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo","histo");
      JetImbalanceJetHemisphere[i][j][0]->hPt = readTH1(RootPath + "hPt_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo","histo");
      JetImbalancePhotonHemisphere[i][j][0]->hPt = readTH1(RootPath+"hPt_imbalance_photon_in_"+(long)(i+1)+ "_Pt_bin_" + (long)(j+1) + "_eta_bin" +DataType + ".root", "histo","histo");
      JetIntrinsic[i][j][0]->hPt = readTH1(RootPath + "hPt_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo","histo");
      
      for(int k=0; k<nAlphaBins; k++){

	JetResponseJetHemisphere[i][j][k]->hResponse = readTH1(RootPath + "response_jet_in_"+ (long)(i+1) +"_Pt_bin_"+ (long)(j+1) +"_eta_bin_"+ (long)(k+1) +"_alpha_bin"+  DataType + ".root","histo","histo");
	JetResponsePhotonHemisphere[i][j][k]->hResponse = readTH1(RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetImbalanceJetHemisphere[i][j][k]->hResponse = readTH1(RootPath + "response_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetImbalancePhotonHemisphere[i][j][k]->hResponse = readTH1(RootPath + "response_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetIntrinsic[i][j][k]->hResponse = readTH1(RootPath + "response_intrinsic_in_" + (long)(i+1) +"_Pt_bin_"+ (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType +".root" , "histo","histo");
	
	
	JetResponseJetHemisphere[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetResponsePhotonHemisphere[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetImbalanceJetHemisphere[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_imbalance_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetImbalancePhotonHemisphere[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_imbalance_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetIntrinsic[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_intrinsic_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	

	// Scale response histograms to the number of entries in MC
	if(JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries() != 0 && JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries() != 0 ){
	  JetResponseJetHemisphere[i][j][k]->hResponse->Scale((JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponseJetHemisphere[i][j][k]->hResponse->Integral()));
	  JetResponsePhotonHemisphere[i][j][k]->hResponse->Scale((JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponsePhotonHemisphere[i][j][k]->hResponse->Integral()));
	  JetImbalancePhotonHemisphere[i][j][k]->hResponse->Scale((JetImbalancePhotonHemisphere[i][j][k]->hResponse->GetEntries())/(JetImbalancePhotonHemisphere[i][j][k]->hResponse->Integral()));
	  JetImbalanceJetHemisphere[i][j][k]->hResponse->Scale((JetImbalanceJetHemisphere[i][j][k]->hResponse->GetEntries())/(JetImbalanceJetHemisphere[i][j][k]->hResponse->Integral()));
	  JetIntrinsic[i][j][k]->hResponse->Scale((JetIntrinsic[i][j][k]->hResponse->GetEntries())/(JetIntrinsic[i][j][k]->hResponse->Integral())); 
	}		
      }
    }
  }
  
  // Remove old object from root_files folder and pdf folder
  TString command;
  command = ".! rm -r " + RootPath + "jet_energy_*" + MethodType + "*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + RootPath + "Resolution_*" + MethodType + "*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + RootPath + "Scale_*" + MethodType + "*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + PDFPath + "jet_energy_*" + MethodType + "*";
  gROOT->ProcessLine(command); 
  command = ".! rm -r " + PDFPath + "Resolution_*" + MethodType + "*";
  gROOT->ProcessLine(command);
  command = ".! rm -r " + PDFPath + "Scale_*" + MethodType + "*";
  gROOT->ProcessLine(command);


  std::cout << "Calculation of Jet energy scale and resolution" << std::endl;

  for(int j=0; j<nEtaBins; j++){
    for(int i=0; i<nPtBins; i++){
      
      JetScaleResAlpha[i][j]  = new CScaleResAlpha;
      JetIntrinsicAlpha[i][j] = new CScaleResAlpha;
      JetImbalanceAlpha[i][j] = new CScaleResAlpha;      
    }
    JetScaleRes[j]         = new CScaleRes;
    JetIntrinsicPt[j]      = new CScaleRes;
    JetImbalancePt[j]      = new CScaleRes;
    
  }
  
  // Fit a Gaussian to all Response Functions or Calculate the RMS (conducted in function calculate()) 
  
  for(int i=0; i<nPtBins; i++){
    for(int j=0; j<nEtaBins; j++){
      
      
      JetResponseJetHemisphere[i][j][0]          -> calculatePt(); 
      JetImbalanceJetHemisphere[i][j][0]         -> calculatePt();
      JetResponsePhotonHemisphere[i][j][0]       -> calculatePt(); 
      JetImbalancePhotonHemisphere[i][j][0]      -> calculatePt();
      JetIntrinsic[i][j][0]                      -> calculatePt(); 
      
      for(int k=0; k<nAlphaBins; k++){

	// Fit Gaussian Functions to Response Histogram if number of entries is large enough 
	if( JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries()>100  && JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries()>100){
	  
	  JetResponseJetHemisphere[i][j][k]       -> calculate(); 
	  JetResponsePhotonHemisphere[i][j][k]    -> calculate();
	  JetImbalanceJetHemisphere[i][j][k]      -> calculate();
	  JetImbalancePhotonHemisphere[i][j][k]   -> calculate();
	  JetIntrinsic[i][j][k]                   -> calculate(); 

	  hChiSquareResolution -> Fill(JetResponseJetHemisphere[i][j][k] -> chi2_ndof,1);
	  hChiSquareResolution -> Fill(JetResponsePhotonHemisphere[i][j][k] -> chi2_ndof,1);
	  hChiSquareImbalance  -> Fill(JetImbalanceJetHemisphere[i][j][k] -> chi2_ndof,1);
	  hChiSquareImbalance  -> Fill(JetImbalancePhotonHemisphere[i][j][k] -> chi2_ndof,1);
	  hChiSquareIntrinsic  -> Fill(JetIntrinsic[i][j][k]->chi2_ndof,1);
	  
        }
	
      }
    }
  }

  extrapolate();

  plotControlPlotsInCalcScale();
 
}





