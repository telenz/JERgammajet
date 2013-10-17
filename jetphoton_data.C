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
  if(date == 2012 && numTrigger!=7){
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

  for(int j=0; j<nEtaBins; j++){
    for(int i=0; i<nPtBins; i++){
      for(int k=0; k<nAlphaBins; k++){
	
	JetResponseJetHemisphere[i][j][k]          = new CResponse(1); 	
	JetResponsePhotonHemisphere[i][j][k]       = new CResponse(1); 	      
      }   
    }
  }
  //-----------------------------------------------------------------------------
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
    for(int k=0; k<nAlphaBins; k++){
      if(alpha >= alphaBins[k] && alpha < alphaBins[k+1]){ 
	
	for(int j=0; j<nEtaBins; j++){
	  
	  if(std::abs(jetEta[corrJets.idx(idx1stJet)])>= etaBins[j] && std::abs(jetEta[corrJets.idx(idx1stJet)]) < etaBins[j+1] ){
	    
	    for(int i=0; i<nPtBins; i++){ 
              
              if(photonPt[0] >= ptBins[i] && photonPt[0] < ptBins[i+1]){

		float deltaphi2ndJet1stJet = 0;
		float deltaphi2ndJetPhoton = 0;
		if(idx2ndJet != -1){
		 
		  deltaphi2ndJet1stJet = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(idx2ndJet)]-jetPhi[corrJets.idx(idx1stJet)]));
		  deltaphi2ndJetPhoton = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(idx2ndJet)]-photonPhi[0])); 
		
	       	}
						 
		if(deltaphi2ndJetPhoton > deltaphi2ndJet1stJet){
		  JetResponseJetHemisphere[i][j][k] -> hResponse      -> Fill(response,weight*PUWeight);
		  JetResponseJetHemisphere[i][j][0] -> hPt            -> Fill(photonPt[0],weight*PUWeight);
		  JetResponseJetHemisphere[i][j][k] -> hAlpha         -> Fill(alpha,weight*PUWeight);
		}
		else{
		  JetResponsePhotonHemisphere[i][j][k] -> hResponse -> Fill(response,weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][0] -> hPt       -> Fill(photonPt[0],weight*PUWeight);
		  JetResponsePhotonHemisphere[i][j][k] -> hAlpha    -> Fill(alpha,weight*PUWeight);
		}
		cutflow<<"nocut = "<<nocut<<endl;
		break;
	      }
	    }
	    break;
	  }
	}
	break;
      }
    }
    
    
    // Fill 2d histogram with first and second Photon of one event
    hPhoton12Pt->Fill(photonPt[0],photonPt[1],weight);
    
    nocut = nocut +1;

  } // End of loop over entries


  TString TotFilename;

  for(int i=0;i<nPtBins;i++){
    for(int j=0; j<nEtaBins; j++){

      // Photon Pt Histograms
      saveObject(JetResponseJetHemisphere[i][j][0]->hPt, RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo");   
      saveObject(JetResponsePhotonHemisphere[i][j][0]->hPt, RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo"); 
      

      for(int k=0; k<nAlphaBins; k++){	

	// Alpha Histograms
	saveObject(JetResponseJetHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_jet_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root","histo");
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hAlpha,RootPath+"hAlpha_photon_in_"+(long)(i+1)+"_Pt_bin_"+(long)(j+1)+"_eta_bin_"+(long)(k+1)+"_alpha_bin"+DataType+".root","histo");
	//Response Histograms
	saveObject(JetResponseJetHemisphere[i][j][k]->hResponse, RootPath + "response_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	saveObject(JetResponsePhotonHemisphere[i][j][k]->hResponse, RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo");
	
      }
    }
  }

  // Save some control plots
  plotControlPlots();


  // CUTFLOW is written into a file
  cutflow.open(RootPath + "cutflow.txt");

  cout<<"CUTFLOW (Data):"<<endl<<"cut1 = "<<cut1<<endl<<"cut2 = "<<cut2<<endl<<"cut3 = "<<cut3<<endl<<"cut4 = "<<cut4<<endl<<"cut5 = "<<cut5<<endl<<"cut6 = "<<cut6<<endl<<"cut7 = "<<cut7<<endl<<"cut8 = "<<cut8<<endl<<"cut9 = "<<cut9<<endl<<"cut10 = "<<cut10<<endl<<"cut11 = "<<cut11<<endl;
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
 

  cout<<"no2ndJetinEvent = "<<no2ndJetinEvent<<endl<<endl;
  cutflow.close();
  delete chain;
  
  filestr.close();   
  //EventVariables.close();
  
}


// --------------------------------------------------------------------------------


// Calculate sigma and mean for all Response functions (also for different nVtx regions)

void calcScale(){


  std::cout << "Read Response Histograms!" << std::endl<< std::endl;

  initializeControlPlotsInCalcScale();

  for(int j=0; j<nEtaBins; j++){
    for(int i=0; i<nPtBins; i++){
      for(int k=0; k<nAlphaBins; k++){

	JetResponseJetHemisphere[i][j][k]         = new CResponse(1); 		
	JetResponsePhotonHemisphere[i][j][k]      = new CResponse(1); 	
      }
    }
  }


  TString TotFilename;  
  for(int i=0;i<nPtBins;i++){
    for(int j=0; j<nEtaBins; j++){
      
      JetResponseJetHemisphere[i][j][0]->hPt = readTH1(RootPath + "hPt_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo","histo");
      JetResponsePhotonHemisphere[i][j][0]->hPt = readTH1(RootPath + "hPt_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin" +  DataType + ".root", "histo","histo");

      for(int k=0; k<nAlphaBins; k++){
	
	JetResponseJetHemisphere[i][j][k]->hResponse = readTH1(RootPath + "response_jet_in_"+ (long)(i+1) +"_Pt_bin_"+ (long)(j+1) +"_eta_bin_"+ (long)(k+1) +"_alpha_bin"+  DataType + ".root","histo","histo");
	JetResponsePhotonHemisphere[i][j][k]->hResponse = readTH1(RootPath + "response_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	
	
	if(testClosure || QCDUncertaintyEvaluation){
	  //Scale the histograms
	  JetResponseJetHemisphere[i][j][k]->hResponse->Scale((JetResponseJetHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponseJetHemisphere[i][j][k]->hResponse->Integral()));
	  JetResponsePhotonHemisphere[i][j][k]->hResponse->Scale((JetResponsePhotonHemisphere[i][j][k]->hResponse->GetEntries())/(JetResponsePhotonHemisphere[i][j][k]->hResponse->Integral()));
	}
	
	JetResponseJetHemisphere[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_jet_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
	JetResponsePhotonHemisphere[i][j][k]->hAlpha = readTH1(RootPath + "hAlpha_photon_in_" + (long)(i+1) + "_Pt_bin_" + (long)(j+1) + "_eta_bin_" + (long)(k+1) + "_alpha_bin" +  DataType + ".root", "histo","histo");
				
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
  

  for(int j=0; j<nEtaBins; j++){
    for(int i=0; i<nPtBins; i++){
      
      JetScaleResAlpha[i][j] = new CScaleResAlpha;
    }
    JetScaleRes[j] = new CScaleRes;
  }
  
  
  std::cout << "Calculation of Jet energy scale and resolution" << std::endl;
    
  // Fit a Gaussian to all Response Functions (conducted in function calculate()) 
  for(int i=0; i<nPtBins; i++){
    for(int j=0; j<nEtaBins; j++){
      
      JetResponseJetHemisphere[i][j][0]          -> calculatePt();
      JetResponsePhotonHemisphere[i][j][0]       -> calculatePt();
      
      for(int k=0; k<nAlphaBins; k++){

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

