// ----------------------------------------------------------------------------------------------
//  Script to calculate JET Energy Response in Gamma Jet Events
//
//  For questions on ROOT see also the manual at
//  http://root.cern.ch/drupal/content/users-guide
// ----------------------------------------------------------------------------------------------

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
#include "TLatex.h"


const bool isMC    = false;

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


// --------- Function and Class declarations -----------------
// -----------------------------------------------------------
void calcScale();
void calcSample();


// Run main script
// This is the actual command to be typed
// used in the ROOT session.
// ---------------------------------------------
int jetphoton_data(int nEvts, int step) {
  
  time_t start,end;
  double dif;
  time (&start);

  // -------------------------------------------------------------------------------------
  if(date == 2011 && numTrigger!=5){
    cout<<"Please correct the number of Triggers used (in CODE/myDeclarations.h)!"<<endl;  
    return 0;
  }
  if(date == 2012 && numTrigger!=8){
    cout<<"Please correct the number of Triggers used (in CODE/myDeclarations.h)!"<<endl;  
    return 0;
  }
  //-------- print out of Information ----------------------------------------------------
  cout<<endl;
  
  if(!isMC)          cout<<"Data!!"<<endl<<endl;
  else if(isMC)      cout<<"MC!!"<<endl<<endl;
  
  if(date == 2012)       cout<<"date = 2012"<<endl<<endl;
  else if(date == 2011)  cout<<"date = 2011"<<endl<<endl;
  
  if(jetType == 1)      cout<<"jetType = ak5PF Jets "<<endl<<endl<<endl;
  else if(jetType == 2) cout<<"jetType = ak5PFCHS Jets"<<endl<<endl<<endl;
  else if(jetType == 4) cout<<"jetType = ak7PFCHS Jets"<<endl<<endl<<endl;
  
  cout<<endl;
  //--------------------------------------------------------------------------------------
  cout<<"calcPath() is executed!"<<endl<<endl;
  calcPath(step);
  cout<<"bookHistos() is executed!"<<endl<<endl;
  bookHistos();
  cout<<"readGammaJet() is executed!"<<endl<<endl;
  readGammaJet(nEvts);
  if(step == 1){
    cout<<"calcSample() is executed!"<<endl<<endl;
    calcSample();
  }
  if(step == 2){
    cout<<"calcScale() is executed!"<<endl<<endl;

    if(detJER == 1)      cout<<endl<<"Method is : Gaus 2 sigma!!"<<endl;
    else if(detJER == 2) cout<<endl<<"Method is : RMS 95!!"<<endl;
    else if(detJER == 3) cout<<endl<<"Method is : RMS 99!!"<<endl;

    calcScale();
  }
  //cout<<"draw() is executed!"<<endl<<endl;
  //draw();
  //--------------------------------------------------------------------------------------

  time (&end);
  dif = difftime(end,start);
  cout<<"elapsed time:"<<dif/60<<" minutes"<<endl<<endl;

  return 0;
}


// Read GammaJet from file and fill histograms. 
// -------------------------------------------------------------------------------
void calcSample() {

  initializeControlPlots();

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      for(int k=0; k<alpha_int; k++){
	
	JetResponseJetHemisphere[i][j][k]          = new CResponse(1); 	
	JetResponsePhotonHemisphere[i][j][k]       = new CResponse(1); 	      
      }   
    }
  }
  //-----------------------------------------------------------------------------

  /*
  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){ 
      h2ndgenJetPtle12GeV[i][j] = createTH1(h2ndgenJetPtle12GeV[i][j],""," 2nd gen Jet pt for jets < 12 GeV ",1000,0,12,"genJetPt");
      TString Filename = HistogramsSmall2ndJetPt + "h2ndgenJetPtle12GeV_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin.root";
      h2ndgenJetPtle12GeV[i][j] = readTH1(Filename, "histo", "histo");
    }
  }
  */
  // Loop over nMax entries, read variables and fill histogram
  
  std::cout << "Processing events" << std::endl;
  int testCalc = 0;
  
  for(int n = 0; n < nMax; n++) {
    
    if( (n+1)%1000000 == 0 ) std::cout << "Event " << (n+1) << std::endl;
    
    // Get this event i.e. entry n and fill values into variables
    chain->GetEntry(n);
    
    //Just MC has PUweights; in data you have to set them to One
    PUWeight= 1.;
    //weight  = 1.;  
  
    //--------------------------------------------------------------------------
    // Sort Jets and apply Cuts
    bool testVar = applyCuts();
    if(!testVar){
      testCalc += 1;
      continue;
    }
    //--------------------------------------------------------------------------
    
    
    // Fill Response Functions    
    for(int k=0; k<alpha_int; k++){
      if(alpha >= alphaBin[k] && alpha < alphaBin[k+1]){ 
	
	for(int j=0; j<eta_int; j++){
	  
	  if(std::abs(jetEta[corrJets.idx(lead_jet)])>= etaBin[j] && std::abs(jetEta[corrJets.idx(lead_jet)]) < etaBin[j+1] ){
	    
	    for(int i=0; i<pt_int; i++){ 
              
              if(photonPt[0] >= bd[i] && photonPt[0] < bd[i+1]){

		
		float deltaphi1stJetPhoton = 0;
		float deltaphi2ndJetPhoton = 0;
		if(jet_2 != -1){
		  if(std::abs(TVector2::Phi_mpi_pi((jetPhi[corrJets.idx(lead_jet)]+photonPhi[0])/2. - jetPhi[corrJets.idx(jet_2)])) < TMath::Pi()/2.){
		    
		    deltaphi1stJetPhoton = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(lead_jet)]-photonPhi[0]));
		  }
		  else deltaphi1stJetPhoton = TMath::Pi()*2. - std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(lead_jet)]-photonPhi[0]));
		  
		  deltaphi2ndJetPhoton = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(jet_2)]-photonPhi[0])); 
		}
		
				 
		if(deltaphi2ndJetPhoton > deltaphi1stJetPhoton/2. || jet_2 == -1){
		  JetResponseJetHemisphere[i][j][k] -> hResponse      -> Fill(response,weight*PUWeight);
		  JetResponseJetHemisphere[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);
		  JetResponseJetHemisphere[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		}
		else{
		  JetResponsePhotonHemisphere[i][j][k] -> hResponse -> Fill(response,weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);
		}
		
		break;
	      }
	      else{
		//cout<<"i = "<<i<<endl;              
                //cout<<"photonPt = "<<photonPt[0]<<endl;
		//cout<<"Something wrong with Trigger and Pt Bin bounds (MC)"<<endl<<endl;
	      }
	    }
	    break;
	  }
	}
	break;
      }
    }
    
    
    // SMALL RESPONSE DIAGNOSTIC::Plots just for small response region
    //if( corrJets.pt(lead_jet)/photonPt[0] <= 0.3 ){
    //hPhotonEta_low_resp_reg -> Fill(photonEta[0],weight);
    
    //if(!filestr.is_open()) cout<<"wrong"<<endl;
    //filestr<<runNum<<":"<<lumiNum<<":"<<eventNum<<endl;	   
    //if(!EventVariables.is_open()) cout<<"EventVariable.txt is not open"<<endl;
    //EventVariables<<photonPt[0]<<":"<<corrJets.pt(lead_jet)<<":"<<photonEta[0]<<":"<<jetEta[corrJets.idx(lead_jet)]<<":"<<photonPhi[0]<<":"<<jetPhi[corrJets.idx(lead_jet)]<<":"<<photonIsoEcal[0]<<":"<<photonIsoHcal[0]<<":"<<photonIsoTrk[0]<<endl;
    //}
    
    //if(( hltPhoton[5] ||  hltPhoton[6] || hltPhoton[7])  && corrJets.pt(lead_jet)/photonPt[0]<=1.0 && photonPt[0]>130.){
    //  hSmallJetPtResponse -> Fill(corrJets.pt(lead_jet),corrJets.pt(lead_jet)/photonPt[0],weight);
    //  hPhotonEta_high_pt_reg_Response -> Fill(photonEta[0],corrJets.pt(lead_jet)/photonPt[0],weight);      
    //}
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
    
    /*// Fill 2d Diagram Response L1 correctio
      hResponseL1 -> Fill(corrJets.pt(lead_jet)/photonPt[0],jetCorrL1[corrJets.idx(lead_jet)],weight);
      hResponseL2L3 -> Fill(corrJets.pt(lead_jet)/photonPt[0],jetCorrL2L3[corrJets.idx(lead_jet)],weight);
      
      //Fill 2d histogram Response agaist Photon Eta
      hRespPhotonEta -> Fill(corrJets.pt(lead_jet)/photonPt[0],photonEta[0],weight);
      
      // Fill 2d Diagram Response against JetEMF JetFHPD JetFRBX
      hResponseJetEMF -> Fill(response,jetEMF[corrJets.idx(lead_jet)],weight);
      hResponseJetFHPD -> Fill(response,jetFHPD[corrJets.idx(lead_jet)],weight);
      hResponseJetFRBX -> Fill(response,jetFRBX[corrJets.idx(lead_jet)],weight);
      
      // Fill 2d Diagram Response against PhotonEMF PhotonFHPD PhotonFRBX
      if(photonidx!=-1){
      hResponsePhotonEMF -> Fill(response,jetEMF[corrJets.idx(photonidx)],weight);
      hResponsePhotonFHPD -> Fill(response,jetFHPD[corrJets.idx(photonidx)],weight);
      hRespoxnsePhotonFRBX -> Fill(response,jetFRBX[corrJets.idx(photonidx)],weight);
      if(response<=1 && corrJets.pt(photonidx)/photonPt[0]<=2)    hResponsePhotonPtRatio -> Fill(response,corrJets.pt(photonidx)/photonPt[0],weight);
      }*/
    
    // Fill 2d histogram with first and second Photon of one event
    hPhoton12Pt->Fill(photonPt[0],photonPt[1],weight);
    
    nocut = nocut +1;

  } // End of loop over entries


  TString TotFilename;

  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){


      saveObject(JetResponseJetHemisphere[i][j][0]->hPt, RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo");   
      saveObject(JetResponsePhotonHemisphere[i][j][0]->hPt, RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo"); 

      for(int k=0; k<alpha_int; k++){	

	saveObject(JetResponseJetHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_jet_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root","histo");
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_photon_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root","histo");
	saveObject(JetResponseJetHemisphere[i][j][k]->hResponse, RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hResponse, RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	
      }
    }
  }

  // Save some control plots
  plotControlPlots();


  // CUTFLOW is written into a file
  cutflow.open(RootPath + "cutflow.txt");

  cout<<"CUTFLOW:"<<endl<<"cut1 = "<<cut1<<endl<<"cut2 = "<<cut2<<endl<<"cut3 = "<<cut3<<endl<<"cut4 = "<<cut4<<endl<<"cut5 = "<<cut5<<endl<<"cut6 = "<<cut6<<endl<<"cut7 = "<<cut7<<endl<<"cut8 = "<<cut8<<endl<<"cut9 = "<<cut9<<endl<<"cut10 = "<<cut10<<endl<<"cut11 = "<<cut11<<endl;
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

  cout<<"number of all events = "<<cut1+cut2+cut3+cut4+cut5+cut6+cut7+cut8+cut9+cut10+cut11+cut12[0]+cut12[1]+cut12[2]+cut12[3]+cut12[4]+cut12[5]+cut12[6]+cut12[7]+cut13[0]+cut13[1]+cut13[2]+cut14+cut15+nocut<<endl<<endl;
  cutflow<<"number of all events = "<<cut1+cut2+cut3+cut4+cut5+cut6+cut7+cut8+cut9+cut10+cut11+cut12[0]+cut12[1]+cut12[2]+cut12[3]+cut12[4]+cut12[5]+cut12[6]+cut12[7]+cut13[0]+cut13[1]+cut13[2]+cut14+cut15+nocut<<endl<<endl;
 


  cout<<endl<<"count2nd = "<<count2nd<<endl<<endl;
  cout<<endl<<"all2ndJets = "<<all2ndJets<<endl<<endl;
  cutflow.close();

  delete chain;
  
  filestr.close();   
  //EventVariables.close();

 
  
}


// --------------------------------------------------------------------------------
// Now take the response functions and fit to the core region a gauss function 
// so that one gets the jet energy scale and jet energy resolution


// Calculate sigma and mean for all Response functions (also for different nVtx regions)

void calcScale(){


  std::cout << "Read Response Histograms!" << std::endl<< std::endl;

  initializeControlPlotsInCalcScale();

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      for(int k=0; k<alpha_int; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 		
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 	
      }
    }
  }


  TFile *file;
  TString TotFilename;  
  for(int i=0;i<pt_int;i++){
    for(int j=0; j<eta_int; j++){

        
      TotFilename = RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";    
      file = new TFile(TotFilename);      
      JetResponseJetHemisphere[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetResponseJetHemisphere[i][j][0]->hPt -> SetDirectory(0);
      delete file;  
      
      TotFilename = RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root";    
      file = new TFile(TotFilename);      
      JetResponsePhotonHemisphere[i][j][0]->hPt =  (TH1D*) gDirectory->Get("histo");
      JetResponsePhotonHemisphere[i][j][0]->hPt -> SetDirectory(0);
      delete file;       

      for(int k=0; k<alpha_int; k++){
	
	TotFilename = RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponseJetHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetResponseJetHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file; 
				
	TotFilename = RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponsePhotonHemisphere[i][j][k]->hResponse =  (TH1D*) gDirectory->Get("histo");
	JetResponsePhotonHemisphere[i][j][k]->hResponse -> SetDirectory(0);
	delete file; 

	if(testClosure){
	  //Scale the histograms
	  JetResponseJetHemisphere[i][j][k]->hResponse->Scale((JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponseJetHemisphere[i][j][k]->hResponse->Integral()));
	  JetResponsePhotonHemisphere[i][j][k]->hResponse->Scale((JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponsePhotonHemisphere[i][j][k]->hResponse->Integral()));
	}
	  		
	TotFilename = RootPath + "hAlpha_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponseJetHemisphere[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetResponseJetHemisphere[i][j][k]->hAlpha -> SetDirectory(0);
	delete file; 
     
	TotFilename = RootPath + "hAlpha_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root";    
	file = new TFile(TotFilename);      
	JetResponsePhotonHemisphere[i][j][k]->hAlpha =  (TH1D*) gDirectory->Get("histo");
	JetResponsePhotonHemisphere[i][j][k]->hAlpha -> SetDirectory(0);
	delete file; 
		
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
  

  for(int j=0; j<eta_int; j++){
    for(int i=0; i<pt_int; i++){
      
      JetScaleResAlpha[i][j] = new CScaleResAlpha;
    }
    JetScaleRes[j] = new CScaleRes;
  }
  
  
  std::cout << "Calculation of Jet energy scale and resolution" << std::endl;
    
  // Fit a Gaussian to all Response Functions (conducted in function calculate()) 
  for(int i=0; i<pt_int; i++){
    for(int j=0; j<eta_int; j++){

      JetResponseJetHemisphere[i][j][0]          -> calculatePt();
      JetResponsePhotonHemisphere[i][j][0]       -> calculatePt();
      
      for(int k=0; k<alpha_int; k++){

	//cout<<"JetResponseJetHemisphere["<<i<<"]["<<j<<"]["<<k<<"]->hResponse->GetEntries() = "<<JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries()<<endl;
	//cout<<"JetResponsePhotonHemisphere["<<i<<"]["<<j<<"]["<<k<<"]->hResponse->GetEntries() = "<<JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries()<<endl;
	
	if(JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries()>100  && JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries()>100){
	  
	  JetResponseJetHemisphere[i][j][k]       -> calculate(); 
	  JetResponsePhotonHemisphere[i][j][k]    -> calculate();

	  hChiSquareResolution -> Fill(JetResponseJetHemisphere[i][j][k] -> chi2_ndof,1);
	  hChiSquareResolution -> Fill(JetResponsePhotonHemisphere[i][j][k] -> chi2_ndof,1);
	  
        }  
	
      }
    }
  }


  extrapolate();  
  plotControlPlotsInCalcScale();

  
}

