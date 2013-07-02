#include "myDeclarations.h"
#include "myClasses.h"
#include "TVector2.h"
#include "TMath.h"


void extrapolate(){

  int idx,length;
  TString etaRegion[nEtaBins], ptRegion[nPtBins], filename, graphTitle, tot_filename, ResImbFilenameScale, ResImbFilenameResolution;
  TFile *f;
  TGraphErrors* ResImbScale[nEtaBins] = {0};
  TGraphErrors* ResImbRes[nEtaBins]   = {0};  
  double_t *YResImbRes      = 0;
  double_t *XResImbRes      = 0;
  double_t *YResImbScale    = 0;
  double_t *YResImbResError = 0; 
  
  // Some description declarations
  for(int j=0;j<nEtaBins;j++){
    etaRegion[j].Form("%4.1f < |#eta| < %4.1f ",etaBins[j],etaBins[j+1]);  
  }

  for(int i=0;i<nPtBins;i++){
    if(i!=nPtBins-1)  ptRegion[i].Form("and %4.0f GeV < p_{T}^{#gamma} < %4.0f GeV",ptBins[i],ptBins[i+1]);
    else             ptRegion[i].Form("and %4.0f GeV < p_{T}^{#gamma} ",ptBins[i]);
  }


  if(!isMC){

    for(int j=0; j<nEtaBins; j++){
      if(date==2012){
      
	if(jetType == 1){
	  ResImbFilenameScale       = "plots_2012/PF_L1FastJet/mc/root_files/Residual_imbalance_Scale_for_";
	  ResImbFilenameResolution  = "plots_2012/PF_L1FastJet/mc/root_files/Residual_imbalance_Resolution_for_";
	}
	else if(jetType == 2  || jetType == 4){
	  ResImbFilenameScale       = "plots_2012/PF_L1CHS/mc/root_files/Residual_imbalance_Scale_for_";
	  ResImbFilenameResolution  = "plots_2012/PF_L1CHS/mc/root_files/Residual_imbalance_Resolution_for_";
	}
      
      }
      else if(date==2011){
	ResImbFilenameScale  = "plots_2011/PF_L1FastJet/mc/root_files/Residual_imbalance_Scale_for_";
	ResImbFilenameResolution  = "plots_2011/PF_L1FastJet/mc/root_files/Residual_imbalance_Resolution_for_";
      }
      ResImbFilenameScale += j+1;
      ResImbFilenameResolution += j+1;

      if(jetType == 1){
	ResImbFilenameScale += "_eta_bin_all_pt_bins_PF_mc";
	ResImbFilenameResolution += "_eta_bin_all_pt_bins_PF_mc";
      }
      else if(jetType == 2 || jetType == 4){
	ResImbFilenameScale += "_eta_bin_all_pt_bins_PFCHS_mc";
	ResImbFilenameResolution += "_eta_bin_all_pt_bins_PFCHS_mc";
      }

      ResImbFilenameScale      += MethodType + ".root";
      ResImbFilenameResolution += MethodType + ".root";

      ResImbScale[j]  = readTGraphErrors(ResImbFilenameScale,"Graph;1","Graph;1");  
      ResImbRes[j]    = readTGraphErrors(ResImbFilenameResolution,"Graph;1","Graph;1");  

      plotTGraphErrors(ResImbRes[j], Form("Residual_imbalances_for_%i_eta_bin",j+1), Form("Residual_imbalances_for_%i_eta_bin",j+1),"p_{T}{#gamma} [GeV]","Residual Imbalance",0.0, 600.,0.,0.5,"",0);
      

      
    }
    
    cout<<"Residual imbalances are read from following file: "<<ResImbFilenameResolution<<endl<<endl;
  }

  // Write in other object all relevant numbers (only if they are different from zero)
  for(int j=0; j<nEtaBins; j++){

    if(!isMC){
      YResImbScale    = ResImbScale[j] -> GetY();        
      YResImbRes      = ResImbRes[j]   -> GetY();
      XResImbRes      = ResImbRes[j]   -> GetX();
      YResImbResError = ResImbRes[j]   -> GetEY();
    }

    for(int i=0; i<nPtBins; i++){
  
      length = nAlphaBins;     
      idx = 0;      
      for(int k=0; k<nAlphaBins;k++){
	
        if(JetResponseJetHemisphere[i][j][k]->mean_array == 0 ){
	  length = length - 1;  
	  continue;
	}
	// Calculate the resulting Resolution from the two histograms in photon and jet hemissphere and write it into 
	// Calculate weighted mean  
	float weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> sigma_error,2);
	//weightPhoton = 0;
	float weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> sigma_error,2);	
	//weightJet = 0;
	float sigma = ((JetResponseJetHemisphere[i][j][k]->sigma_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->sigma_array)*weightPhoton)/(weightJet+weightPhoton);
	float sigmaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][k] -> alpha_error,2);
	//weightPhoton = 0;
	weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][k] -> alpha_error,2);	
	//weightJet = 0 ;
	float alphaVal = ((JetResponseJetHemisphere[i][j][k]->alpha_array)*weightJet + (JetResponsePhotonHemisphere[i][j][k]->alpha_array)*weightPhoton)/(weightJet+weightPhoton);
	float alphaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	
	JetScaleResAlpha[i][j]  -> alpha[idx]      = alphaVal;
        JetScaleResAlpha[i][j]  -> alphaError[idx] = alphaError;
	JetScaleResAlpha[i][j]  -> sigma[idx]      = sigma;         
        JetScaleResAlpha[i][j]  -> sigmaError[idx] = sigmaError;	
	JetScaleResAlpha[i][j]  -> mean[idx]       = JetResponseJetHemisphere[i][j][k] -> mean_array;         
        JetScaleResAlpha[i][j]  -> meanError[idx]  = JetResponseJetHemisphere[i][j][k] -> mean_error;
	

	// Only for MC important (imbalance and intrinsic response)
	if(isMC){
	  weightPhoton = 1/TMath::Power(JetImbalancePhotonHemisphere[i][j][k] -> sigma_error,2);
	  //weightPhoton = 0;
	  weightJet = 1/TMath::Power(JetImbalanceJetHemisphere[i][j][k] -> sigma_error,2);	
	  //weightJet = 0;
	  sigma = ((JetImbalanceJetHemisphere[i][j][k]->sigma_array)*weightJet + (JetImbalancePhotonHemisphere[i][j][k]->sigma_array)*weightPhoton)/(weightJet+weightPhoton);
	  sigmaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	  weightPhoton = 1/TMath::Power(JetImbalancePhotonHemisphere[i][j][k] -> alpha_error,2);
	  //weightPhoton = 0;
	  weightJet = 1/TMath::Power(JetImbalanceJetHemisphere[i][j][k] -> alpha_error,2);
	  //weightJet = 0;	
	  alphaVal = ((JetImbalanceJetHemisphere[i][j][k]->alpha_array)*weightJet + (JetImbalancePhotonHemisphere[i][j][k]->alpha_array)*weightPhoton)/(weightJet+weightPhoton);
	  alphaError = TMath::Sqrt(1./(weightJet+weightPhoton));
	  
	  JetImbalanceAlpha[i][j] -> alpha[idx]      = alphaVal;  
	  JetImbalanceAlpha[i][j] -> alphaError[idx] = alphaError; 
	  JetImbalanceAlpha[i][j] -> sigma[idx]      = sigma;             
	  JetImbalanceAlpha[i][j] -> sigmaError[idx] = sigmaError;   
	  JetImbalanceAlpha[i][j] -> mean[idx]       = JetImbalanceJetHemisphere[i][j][k] -> mean_array;         
	  JetImbalanceAlpha[i][j] -> meanError[idx]  = JetImbalanceJetHemisphere[i][j][k] -> mean_error;

	  
	  JetIntrinsicAlpha[i][j] -> alpha[idx]      = JetIntrinsic[i][j][k] -> alpha_array;  
	  JetIntrinsicAlpha[i][j] -> alphaError[idx] = JetIntrinsic[i][j][k] -> alpha_error;
	  JetIntrinsicAlpha[i][j] -> mean[idx]       = JetIntrinsic[i][j][k] -> mean_array;         
	  JetIntrinsicAlpha[i][j] -> meanError[idx]  = JetIntrinsic[i][j][k] -> mean_error;       
	  JetIntrinsicAlpha[i][j] -> sigma[idx]      = JetIntrinsic[i][j][k] -> sigma_array;         
	  JetIntrinsicAlpha[i][j] -> sigmaError[idx] = JetIntrinsic[i][j][k] -> sigma_error;
	}
	idx = idx + 1;  
	
      }



      // Calculate weighted mean of Photon Pt
      float weightPhoton = 1/TMath::Power(JetResponsePhotonHemisphere[i][j][0] -> pT_error,2);
      //weightPhoton = 0;
      float weightJet = 1/TMath::Power(JetResponseJetHemisphere[i][j][0] -> pT_error,2);	
      float pT = ((JetResponseJetHemisphere[i][j][0]->pT_array)*weightJet + (JetResponsePhotonHemisphere[i][j][0]->pT_array)*weightPhoton)/(weightJet+weightPhoton);
      //weightJet = 0;
      float pTError = TMath::Sqrt(1./(weightJet+weightPhoton));

      JetScaleResAlpha[i][j]  -> pT        = pT;
      JetScaleResAlpha[i][j]  -> pTError   = pTError;

      if(isMC){
	JetImbalanceAlpha[i][j] -> pT        = pT;  
	JetImbalanceAlpha[i][j] -> pTError   = pTError;
	
	JetIntrinsicAlpha[i][j] -> pT         = JetIntrinsic[i][j][0] -> pT_array;  
	JetIntrinsicAlpha[i][j] -> pTError    = JetIntrinsic[i][j][0] -> pT_error;
      }
      
      	    
      cout<<endl<<endl<<"ptBoundLow = "<<ptBins[i]<<endl;
      cout<<"Eta bin = "<<j+1<<endl;

      if(isMC){
	JetIntrinsicAlpha[i][j] -> calculate(length, 2);   
	JetImbalanceAlpha[i][j] -> calculate(length, 3);
   
	JetScaleResAlpha[i][j]  -> qprime      = JetImbalanceAlpha[i][j] -> qprime;
	JetScaleResAlpha[i][j]  -> qprimeError = JetImbalanceAlpha[i][j] -> qprimeError;
	JetScaleResAlpha[i][j]  -> mprime      = JetImbalanceAlpha[i][j] -> mprime;
	JetScaleResAlpha[i][j]  -> cprime      = JetIntrinsicAlpha[i][j] -> cprime;
	JetScaleResAlpha[i][j]  -> cprimeError = JetIntrinsicAlpha[i][j] -> cprimeError;
	JetScaleResAlpha[i][j]  -> q           = JetImbalanceAlpha[i][j] -> q;
	JetScaleResAlpha[i][j]  -> qError      = JetImbalanceAlpha[i][j] -> qError;
      }
      else if(!isMC){

	// fix q and qprime on residual imbalance value (for systematic Uncertainies take more/less 0.5*q
	JetScaleResAlpha[i][j] -> q           = YResImbScale[i];//-0.25*YResImbScale[i];
	JetScaleResAlpha[i][j] -> qprime      = YResImbRes[i];//-0.25*YResImbRes[i];
	JetScaleResAlpha[i][j] -> qprimeError = YResImbResError[i]  ;
	
	cout<<"qprime["<<i<<"] = "<<YResImbRes[i]<<endl;
	cout<<"Xqprime["<<i<<"] = "<<XResImbRes[i]<<endl;
	cout<<"JetScaleResAlpha[i][j]  -> pT ="<<JetScaleResAlpha[i][j]  -> pT <<endl<<endl;
      }

      JetScaleResAlpha[i][j]  -> calculate(length, 4);
            
      char legEntry[100];
      
	if(JetScaleResAlpha[i][j]->scale !=0){

	  filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin",j+1,i+1);
	  graphTitle = "Scale for " + etaRegion[j] + " " + ptRegion[i];
	  
	  plotTGraphErrors(JetScaleResAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha", "JES",0,alphaBins[nAlphaBins],0.8,1.2, "",0);
	  filename = (TString) "jet_energy_scale_for_" + (long) (j+1) + (TString) "_eta_bin_" + (long) (i+1) + (TString) "_pTGamma_bin_intrinsic";
	  graphTitle = "Scale for " + etaRegion[j] + " " + ptRegion[i] + " (intrinsic)";

	  if(isMC){
	    plotTGraphErrors(JetIntrinsicAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha","JES",0,alphaBins[nAlphaBins],0.8,1.2,"",0);
	    filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_imbalance",j+1,i+1);
	    graphTitle = "Scale for " + etaRegion[j] + " " + ptRegion[i] + " (imbalance)";
	    
	    plotTGraphErrors(JetImbalanceAlpha[i][j] -> gJetScaleAlpha, filename, graphTitle, "alpha","JES",0,alphaBins[nAlphaBins],0.8,1.2,"",0);  

	    
	    TF1* totalScale = new TF1("totalScale"," [0]*(1. - [1] - [2] *TMath::Power(x,2))",0,600); 
	    totalScale -> SetParameters(JetIntrinsicAlpha[i][j]->c, JetImbalanceAlpha[i][j]->q,JetImbalanceAlpha[i][j]->m);
	    filename.Form("jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_total",j+1,i+1);
	    tot_filename = RootPath + filename + DataType + MethodType + ".root";
	    f = new TFile(tot_filename,"RECREATE");
	    f -> WriteTObject(totalScale);
	    f ->Close();
	    delete f;
	  }
        
      
	  filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin",j+1,i+1);
	  graphTitle = "Resolution for " + etaRegion[j] + " " + ptRegion[i];
	  sprintf(legEntry,"#splitline{Chi^2 = %4.2f}{ndof = %i}",JetScaleResAlpha[i][j] -> gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetChisquare(),JetScaleResAlpha[i][j] -> gJetResolutionAlpha->GetFunction("fResolutionAlpha") ->GetNDF());
	  plotTGraphErrors(JetScaleResAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,alphaBins[nAlphaBins],0.,0.4,legEntry,1);

	  if(isMC){
	    filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_intrinsic",j+1,i+1);
	    graphTitle = "Resolution for " + etaRegion[j] + " " + ptRegion[i] + " (intrinsic)";
	    sprintf(legEntry,"#splitline{Chi^2 = %4.2f}{ndof = %i}",JetIntrinsicAlpha[i][j]->gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetChisquare(),JetIntrinsicAlpha[i][j] -> gJetResolutionAlpha->GetFunction("fResolutionAlpha") ->GetNDF());
	    plotTGraphErrors(JetIntrinsicAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,alphaBins[nAlphaBins],0.,0.4,legEntry,1);
	    filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_imbalance",j+1,i+1);
	    graphTitle = "Resolution for " + etaRegion[j] + " " + ptRegion[i] + " (imbalance)";
	    sprintf(legEntry,"#splitline{Chi^2 = %4.2f}{ndof = %i}",JetImbalanceAlpha[i][j]->gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetChisquare(),JetImbalanceAlpha[i][j]->gJetResolutionAlpha->GetFunction("fResolutionAlpha")->GetNDF());
	    plotTGraphErrors(JetImbalanceAlpha[i][j] -> gJetResolutionAlpha, filename, graphTitle, "alpha","JER",0,alphaBins[nAlphaBins],0.,0.4,legEntry,1); 
	  
	    TF1* totalResolution = new TF1("totalResolution","TMath::Sqrt(TMath::Power([0],2) + TMath::Power([1],2) + 2*[1]*[2]*x + TMath::Power(([2]*x),2) )",0,600);
	    totalResolution -> SetParameters(JetIntrinsicAlpha[i][j]->cprime, JetImbalanceAlpha[i][j]->qprime,JetImbalanceAlpha[i][j]->mprime);
	    filename.Form("jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_total",j+1,i+1);  
	    tot_filename = RootPath + filename + DataType + MethodType +".root";
	    f = new TFile(tot_filename,"RECREATE");
	    f -> WriteTObject(totalResolution);
	    f ->Close();
	    delete f; 
	  }
	}
    }
  }
   
  
  // DRAW LAST PLOTS
  char legEntry[100];
  double q[nPtBins],qprime[nPtBins],qError[nPtBins],qprimeError[nPtBins]; 
  double pT[nPtBins],pTError[nPtBins];
  
  
  for(int j=0; j<nEtaBins; j++){
    
    length = nPtBins;
    idx = 0;
    
    for(int i=0; i<nPtBins; i++){


      if(isMC){
	q[i]           = JetImbalanceAlpha[i][j] -> q;
	qError[i]      = JetImbalanceAlpha[i][j] -> qError;
	qprime[i]      = JetImbalanceAlpha[i][j] -> qprime;
	qprimeError[i] = JetImbalanceAlpha[i][j] -> qprimeError;
	pT[i]          = JetScaleResAlpha[i][j]  -> pT;
	pTError[i]     = JetScaleResAlpha[i][j]  -> pTError;
      }

      if(JetScaleResAlpha[i][j] -> scale == 0 ){
        length = length - 1;       
        continue;
      }
 
      JetScaleRes[j] -> pT[idx]              = JetScaleResAlpha[i][j] -> pT;
      JetScaleRes[j] -> pTError[idx]         = JetScaleResAlpha[i][j] -> pTError;
      JetScaleRes[j] -> scale[idx]           = JetScaleResAlpha[i][j] -> scale;
      JetScaleRes[j] -> scaleError[idx]      = JetScaleResAlpha[i][j] -> scaleError;
      JetScaleRes[j] -> resolution[idx]      = JetScaleResAlpha[i][j] -> resolution;
      JetScaleRes[j] -> resolutionError[idx] = JetScaleResAlpha[i][j] -> resolutionError;

      if(isMC){
	JetImbalancePt[j] -> pT[idx]           = JetScaleRes[j] -> pT[idx];
	JetImbalancePt[j] -> pTError[idx]      = JetScaleRes[j] -> pTError[idx];     
	JetIntrinsicPt[j] -> pT[idx]           = JetIntrinsicAlpha[i][j] -> pT;
	JetIntrinsicPt[j] -> pTError[idx]      = JetIntrinsicAlpha[i][j] -> pTError;  
	
	JetScaleRes[j] -> q[idx]               = JetImbalanceAlpha[i][j] -> q;
	JetScaleRes[j] -> qError[idx]          = JetImbalanceAlpha[i][j] -> qError;
	JetScaleRes[j] -> qprime[idx]          = JetImbalanceAlpha[i][j] -> qprime;
	JetScaleRes[j] -> qprimeError[idx]     = JetImbalanceAlpha[i][j] -> qprimeError;
	
	JetIntrinsicPt[j] -> scale[idx]           = JetIntrinsicAlpha[i][j] -> scale;
	JetIntrinsicPt[j] -> scaleError[idx]      = JetIntrinsicAlpha[i][j] -> scaleError;
	JetIntrinsicPt[j] -> resolution[idx]      = JetIntrinsicAlpha[i][j] -> resolution;
	JetIntrinsicPt[j] -> resolutionError[idx] = JetIntrinsicAlpha[i][j] -> resolutionError;
	
	JetImbalancePt[j] -> scale[idx]           = JetImbalanceAlpha[i][j] -> scale;
	JetImbalancePt[j] -> scaleError[idx]      = JetImbalanceAlpha[i][j] -> scaleError;
	JetImbalancePt[j] -> resolution[idx]      = JetImbalanceAlpha[i][j] -> resolution;
	JetImbalancePt[j] -> resolutionError[idx] = JetImbalanceAlpha[i][j] -> resolutionError;  
      }

      idx = idx +1;
      
    }
      
    TGraphErrors *gq          = new TGraphErrors(nPtBins, pT , q, pTError, qError);
    TGraphErrors *gqprime     = new TGraphErrors(nPtBins, pT , qprime, pTError, qprimeError);  

    // Calculate Fit variable for Resolution and Scale
    cout<<"Fit is done for JER(photonPt) (Full Resolution):"<<endl;
    JetScaleRes[j]    -> calculate(length);
    if(isMC){
      cout<<"Fit is done for JER(photonPt) (Intrinsic):"<<endl;
      JetIntrinsicPt[j] -> calculate(length);
      cout<<"Fit is done for JER(photonPt) (Imbalance):"<<endl;
      JetImbalancePt[j] -> calculate(length);  
    }
    
    // "Scale" Plots
    TString filenameSubsample;
    filename.Form("Scale_for_%i_eta_bin",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Scale for " + etaRegion[j];
    sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetScaleRes[j] -> scalenumber,JetScaleRes[j] ->scalenumberError);
    plotTGraphErrors(JetScaleRes[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 

    if(isMC){
      filename.Form("Scale_for_%i_eta_bin_intrinsic",j+1);
      filename += filenameSubsample;
      graphTitle = "Jet Energy Scale for " + etaRegion[j] + " (intrinsic)";
      sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetIntrinsicPt[j] -> scalenumber,JetIntrinsicPt[j] ->scalenumberError);
      plotTGraphErrors(JetIntrinsicPt[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 
      
      filename.Form("Scale_for_%i_eta_bin_imbalance",j+1);
      filename += filenameSubsample;
      graphTitle = "Jet Energy Scale for " + etaRegion[j] + " (imbalance)";
      sprintf(legEntry,"scale = %4.3f #pm %4.3f ",JetImbalancePt[j] -> scalenumber,JetImbalancePt[j] ->scalenumberError);
      plotTGraphErrors(JetImbalancePt[j] -> gScale, filename, graphTitle, "p_{T}^{#gamma}","JES",0,600,0.7,1.1,legEntry,1); 

      filename.Form("Residual_imbalance_Scale_for_%i_eta_bin",j+1);
      filename += filenameSubsample;
      graphTitle = "Residual Imbalance for " + etaRegion[j] + " (Scale)";
      plotTGraphErrors(JetScaleRes[j] -> gq, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);

      filename.Form("Residual_imbalance_Scale_for_%i_eta_bin_all_pt_bins",j+1);
      filename += filenameSubsample;
      graphTitle = "Residual Imbalance for " + etaRegion[j] + " (Scale)";
      plotTGraphErrors(gq, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);
    }
 

    // Resolution plots
    filename.Form("Resolution_for_%i_eta_bin",j+1);
    filename += filenameSubsample;
    graphTitle = "Jet Energy Resolution for " + etaRegion[j];
    sprintf(legEntry,"#splitline{N = %4.3f #pm %4.3f, S = %4.3f #pm %4.3f}{m = %4.3f #pm %4.3f, C = %4.3f #pm %4.3f}",JetScaleRes[j]->N,JetScaleRes[j]->NErr,JetScaleRes[j]->S,JetScaleRes[j]->SErr,JetScaleRes[j]->m,JetScaleRes[j]->mErr,JetScaleRes[j]->C,JetScaleRes[j]->CErr);
    plotTGraphErrors(JetScaleRes[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,legEntry,1); 

    if(isMC){
      filename.Form("Resolution_for_%i_eta_bin_intrinsic",j+1);
      filename += filenameSubsample;
      graphTitle = "Jet Energy Resolution for " + etaRegion[j] + "(intrinsic)";
      sprintf(legEntry,"#splitline{N = %4.3f #pm %4.3f, S = %4.3f #pm %4.3f}{m = %4.3f #pm %4.3f, C = %4.3f #pm %4.3f}",JetIntrinsicPt[j]->N,JetIntrinsicPt[j]->NErr,JetIntrinsicPt[j]->S,JetIntrinsicPt[j]->SErr,JetIntrinsicPt[j]->m,JetIntrinsicPt[j]->mErr,JetIntrinsicPt[j]->C,JetIntrinsicPt[j]->CErr);
      plotTGraphErrors(JetIntrinsicPt[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,legEntry,1); 

      filename.Form("Resolution_for_%i_eta_bin_imbalance",j+1);
      filename += filenameSubsample;
      graphTitle = "Jet Energy Resolution for " + etaRegion[j] + "(imbalance)";
      plotTGraphErrors(JetImbalancePt[j]->gResolution, filename, graphTitle, "p_{T}^{#gamma}","JER",0,600,-0.1,0.3,"a",0); 

      filename.Form("Residual_imbalance_Resolution_for_%i_eta_bin",j+1);
      filename += filenameSubsample;
      graphTitle = "Residual Imbalance for " + etaRegion[j]+ " (Resolution)";
      plotTGraphErrors(JetScaleRes[j] -> gqprime, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);

      filename.Form("Residual_imbalance_Resolution_for_%i_eta_bin_all_pt_bins",j+1);
      filename += filenameSubsample;
      graphTitle = "Residual Imbalance for " + etaRegion[j]+ " (Resolution)";
      plotTGraphErrors(gqprime, filename, graphTitle, "p_{T}^{#gamma}","Residual Imbalance",0,600,-0.04,0.1,"a",0);
    }
    
  } // End of Eta loop
}


  
 
  
   

   
