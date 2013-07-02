#include "myDeclarations.h"
#include "myClasses.h"
#include "myFunctions.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"


    
CResponse::CResponse(int typeResponse){
     
  sigma_array = 0.5;
  sigmaAUX = 0.5;
  sigma_error = 0;
  mean_array = 0.0;
  mean_error = 0;
  alpha_array = 0;
  alpha_error = 0;
  pT_array = 0;
  pT_error = 0;
  RMS_array = 0;
  RMS_error = 0;
  chi2 = 0;
  ndof = 0;
  chi2_ndof = 0;
  meanAUX = 1.;
  meanAUX2 = 1.;
  sigma = 0;
  

  if(typeResponse==1)      hResponse  = createTH1(hResponse,"","Response",5000,0,2,"p_{T}^{recoJet}/p_{T}^{#gamma}");
  else if(typeResponse==2) hResponse  = createTH1(hResponse,"","Intrinsic Response",5000,0,2,"p_{T}^{recoJet}/p_{T}^{genJet}");
  else if(typeResponse==3) hResponse  = createTH1(hResponse,"","Imbalance Response",5000,0,2,"p_{T}^{genJet}/p_{T}^{#gamma}");

  hAlpha = createTH1(hAlpha,"","p_{T}^(Jet_{2})/p_{T}^{#gamma}",10000,0,alphaBins[nAlphaBins],"alpha");
  hPt = createTH1(hPt,"","p_{T}^{#gamma} - Distribution",10000,0,5000,"p^{#gamma}_{T}");
  hPtLeadingJet    = createTH1(hPt,"","p_{T}^{1st jet} - Distribution",10000,0,5000,"p^{1st jet}_{T}");
  hgenPtLeadingJet = createTH1(hPt,"","p_{T}^{1st gen. jet} - Distribution",10000,0,5000,"p^{1st gen. jet}_{T}");
  
};
 

CResponse::~CResponse(){
  delete hResponse;
  delete hAlpha;
  delete hPt;   
};

void CResponse::calculatePt(){
  
  pT_array    = hPt    -> GetMean(1);
  pT_error    = hPt    -> GetMeanError(1);
}

void CResponse::calculate(){
   
  alpha_array = hAlpha -> GetMean(1);
  alpha_error = hAlpha -> GetMeanError(1);
  
  RMS_array = hResponse -> GetRMS(1);
  RMS_error = hResponse -> GetRMSError(1);

  meanAUX  = hResponse -> GetMean(1);
  sigmaAUX = hResponse -> GetRMS(1);

  float eps = 1.0;
  int l = 0;
  TF1* f1 = new TF1("gauss","gaus",0,2); 
  TFitResultPtr r;

  if(detJER == 1){

    while(eps > 0.005 && l <10){
      
      sigma_array = sigmaAUX;
    
      r = hResponse -> Fit("gauss","QL","",meanAUX-2.0*sigma_array,meanAUX+2.0*sigma_array);  
      
      sigmaAUX = f1 -> GetParameter(2);
      meanAUX  = f1 -> GetParameter(1);
      
      if (sigma_array>=sigmaAUX) eps = sigma_array/sigmaAUX-1.;
      else                       eps = sigmaAUX/sigma_array -1.;
      l++;    
    }

    //TMatrixDSym cov = r->GetCovarianceMatrix();

    sigma_array = (f1 -> GetParameter(2))/(f1 -> GetParameter(1));
    // Error taken the correlation between the mean and sigma into account
    //sigma_error = sqrt(TMath::Power(1./f1->GetParameter(1),2)*TMath::Power(f1->GetParError(2),2) + TMath::Power(f1->GetParameter(2)/TMath::Power(f1->GetParameter(1),2),2)*TMath::Power(f1->GetParError(1),2) - 2.0*f1->GetParameter(2)/TMath::Power(f1->GetParameter(1),3)*cov[1][2]);
    sigma_error = sqrt(TMath::Power(1./f1->GetParameter(1),2)*TMath::Power(f1->GetParError(2),2) + TMath::Power(f1->GetParameter(2)/TMath::Power(f1->GetParameter(1),2),2)*TMath::Power(f1->GetParError(1),2));


    mean_array  = f1 -> GetParameter(1);
    mean_error  = f1 -> GetParError(1);
    
    chi2 = f1 -> GetChisquare();
    ndof = f1 -> GetNDF();  
    chi2_ndof = chi2/(1.0*ndof);
 
    //delete f1; 
  }
  else {

    double rmsRange  = 0;
    //double gausSigma = 0;
    
    if(detJER == 2){
      rmsRange  = 0.95;
      //gausSigma = 1.960;
    }
    else if(detJER == 3){
      rmsRange = 0.99;  
      //gausSigma = 2.576;
    }
   

    while(eps > 0.005 && l <3){
      
      sigma_array = sigmaAUX;
    
      r = hResponse -> Fit("gauss","QL","",meanAUX-2.0*sigma_array,meanAUX+2.0*sigma_array);  
      
      sigmaAUX = f1 -> GetParameter(2);
      meanAUX  = f1 -> GetParameter(1);
      
      if (sigma_array>=sigmaAUX) eps = sigma_array/sigmaAUX-1.;
      else                       eps = sigmaAUX/sigma_array -1.;
      l++;    
    }
    
    //TMatrixDSym cov = r -> GetCovarianceMatrix();
    int lowBin  = 0;
    int highBin = 0;
    
    for(int i = 0; i<hResponse->GetNbinsX(); i++){
      
      if(hResponse->Integral(hResponse->FindBin(meanAUX)-i, hResponse->FindBin(meanAUX)+i)/hResponse->Integral() >= rmsRange){
	
	lowBin  =  hResponse->FindBin(meanAUX)-i;
	highBin =  hResponse->FindBin(meanAUX)+i;
	break;
      }
    }
      
    hResponse -> GetXaxis() -> SetRange(lowBin,highBin);

    sigma_array = (hResponse->GetRMS(1))/(hResponse->GetMean(1));
    
    
    //    sigma_error = sqrt(TMath::Power(1./hResponse->GetMean(1),2)*TMath::Power(hResponse->GetRMSError(1),2) + TMath::Power(hResponse->GetRMS(1)/TMath::Power(hResponse->GetMean(1),2),2)*TMath::Power(hResponse->GetMeanError(1),2) - 2.0*f1->GetParameter(2)/TMath::Power(f1->GetParameter(1),3)*cov[1][2]);
    sigma_error = sqrt(TMath::Power(1./hResponse->GetMean(1),2)*TMath::Power(hResponse->GetRMSError(1),2) + TMath::Power(hResponse->GetRMS(1)/TMath::Power(hResponse->GetMean(1),2),2)*TMath::Power(hResponse->GetMeanError(1),2));
    
    
    mean_array = hResponse -> GetMean(1);
    mean_error = hResponse -> GetMeanError(1);
    
  }

};

//----
CScaleResAlpha::CScaleResAlpha(){
  
  pT = 0;
  pTError = 0;  
  scale = 0;
  scaleError = 0;
  resolution = 0;
  resolutionError = 0;
  q = 0;
  c = 0;
  m = 0;
  qprime = 0;
  qprimeError = 0;
  cprime = 0;
  cprimeError = 0;
  mprime = 0;
  mprimeError = 0;
  
  for(int i=0; i<2*nAlphaBins; i++){
    alpha[i] = 0;
    alphaError[i] = 0;
    mean[i] = 0;
    meanError[i] = 0;
    sigma[i] = 0;
    sigmaError[i] = 0; 
  }
}


void CScaleResAlpha::calculate(int length, int fit){
  
 
  if(length>=4){
   
    gJetScaleAlpha = new TGraphErrors(length, alpha , mean, alphaError, meanError);
    gJetResolutionAlpha = new TGraphErrors(length, alpha , sigma, alphaError, sigmaError);
    
    if(fit == 2){  //mc intrinsic
      fScaleAlpha = new TF1("fScaleAlpha","[0]",0,alphaBins[nAlphaBins]);
      gJetScaleAlpha -> Fit("fScaleAlpha","QR");
      c = fScaleAlpha -> GetParameter(0);
      scale      = fScaleAlpha -> GetParameter(0);
      scaleError = fScaleAlpha -> GetParError(0);
      
      fResolutionAlpha = new TF1("fResolutionAlpha","[0]",0,alphaBins[nAlphaBins]); 
      gJetResolutionAlpha -> Fit("fResolutionAlpha","QR");
      cprime      = fResolutionAlpha -> GetParameter(0);
      cprimeError = fResolutionAlpha -> GetParError(0);
      resolution      = fResolutionAlpha -> GetParameter(0);
      resolutionError = fResolutionAlpha -> GetParError(0);

      cout<<"Intrinsic = "<<endl;
      cout<<"cprime from intrinsic response      = "<<cprime<<endl;
      cout<<"cprimeError from intrinsic response = "<<resolutionError<<endl<<endl;      
    }
    else if(fit == 3){   //mc imbalance
      
      fScaleAlpha = new TF1("fScaleAlpha","1. - [0] - [1]*TMath::Power(x,2)",0,alphaBins[nAlphaBins]); 
      gJetScaleAlpha -> Fit("fScaleAlpha","QR");
      q = fScaleAlpha -> GetParameter(0);
      qError = fScaleAlpha -> GetParError(0);
      m = fScaleAlpha -> GetParameter(1);
      scale      = 1 - fScaleAlpha -> GetParameter(0);
      scaleError = fScaleAlpha -> GetParError(0);
     
      fResolutionAlpha = new TF1("fResolutionAlpha","[0] + [1]*x",0,alphaBins[nAlphaBins]); 
      gJetResolutionAlpha -> Fit("fResolutionAlpha","QR");
      qprime      = fResolutionAlpha -> GetParameter(0);
      qprimeError = fResolutionAlpha -> GetParError(0);
      

      if(qprime<0){
	cout<<"qprime is negative !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	cout<<"qprime      = "<<qprime<<endl;
	cout<<"qprimeError = "<<qprimeError<<endl<<endl;;
	  if(qprime + qprimeError > 0){
	    cout<<"But with Error compatible to Zero! - qprime Set to zero"<<endl;
	    qprime      = 0;
	    qprimeError = 0;
	    mprime      = fResolutionAlpha -> GetParameter(1);
	    fResolutionAlpha -> FixParameter(0,0);
	    fResolutionAlpha -> SetParameter(1,mprime);

	    gJetResolutionAlpha -> Fit("fResolutionAlpha","QR");
	  }
	  else{
	    cout<<"!!!!!!!!!!!!!Not compatible with zero including ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (Please check that bin)"<<endl;
	    
	    qprime      = 0;
	    qprimeError = 0;
	    mprime      = fResolutionAlpha -> GetParameter(1);
	    fResolutionAlpha -> FixParameter(0,0);
	    fResolutionAlpha -> SetParameter(1,mprime);
	    
	    gJetResolutionAlpha -> Fit("fResolutionAlpha","QR");
	  }
      }

      mprime          = fResolutionAlpha -> GetParameter(1);

      resolution      = qprime;
      resolutionError = qprimeError;
   
      cout<<"qprime from imbalance response      = "<<qprime<<endl;
      cout<<"qprimeError from imbalance response = "<<qprimeError<<endl<<endl;;

    }
    else if(fit == 4){ //data and MC
      
      fScaleAlpha      = new TF1("fScaleAlpha","[0]*(1. - [1] - [2]*TMath::Power(x,2))",0,alphaBins[nAlphaBins]); 
      fScaleAlpha -> FixParameter(1,q);
      gJetScaleAlpha -> Fit("fScaleAlpha","QR");    
      scale      = fScaleAlpha -> GetParameter(0);
      scaleError = fScaleAlpha -> GetParError(0);

      // Resolution       
      fResolutionAlpha = new TF1("fResolutionAlpha","TMath::Sqrt(TMath::Power([0],2) +TMath::Power([1],2) +2*[1]*[2]*x + TMath::Power(([2]*x),2) )",0,alphaBins[nAlphaBins]); 
      fResolutionAlpha    -> FixParameter(1,qprime);        
      gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");

      
      if(fResolutionAlpha -> GetParameter(2) < 0){	

	  cout<<"----------------------------------------------------Fit Paramter mprime < 0 -> Set bounds on it !!!!"<<endl;
	  cout<<"mprime      = "<<fResolutionAlpha -> GetParameter(2)<<endl;
	  cout<<"mprimeError = "<<fResolutionAlpha -> GetParError(2)<<endl;

	  fResolutionAlpha -> SetParLimits(2,0.0,10000000);
	  fResolutionAlpha -> SetParameter(2,0.005);

	  if(isMC){
	    fResolutionAlpha -> SetParameter(0,cprime);
	  }
	  else if(!isMC){
	    fResolutionAlpha -> SetParameter(0,0.05);
	  }
	  gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");	

	  cout<<"mprime after setting bound:"<<endl;
	  cout<<"mprime      = "<<fResolutionAlpha -> GetParameter(2)<<endl;
	  cout<<"mprimeError = "<<fResolutionAlpha -> GetParError(2)<<endl;

	  fResolutionAlpha -> FixParameter(2,fResolutionAlpha -> GetParameter(2));
	  cout<<"mprime is now fixed to that value to avoid strange errors of cmprime!"<<endl;
	  gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");	
	  
	
      }

      if(fResolutionAlpha -> GetParameter(0) < 0){

	cout<<"--------------------------------------------------- Fit Paramter cprime < 0 -> Set bounds on it!!!!"<<endl;

	//fResolutionAlpha -> SetParLimits(2,0,10000000);
	fResolutionAlpha -> SetParLimits(0,0,10000000);
	
	if(isMC){
	  fResolutionAlpha -> SetParameter(0,cprime);
	  fResolutionAlpha -> SetParameter(2,mprime);
	}
	else if(!isMC){
	  fResolutionAlpha -> SetParameter(0,0.05);
	  fResolutionAlpha -> SetParameter(2,0.001);
	}
	
	gJetResolutionAlpha -> Fit("fResolutionAlpha","QRB");	

      }
      

      resolution      = fResolutionAlpha -> GetParameter(0);
      resolutionError = fResolutionAlpha -> GetParError(0);
      
      cout<<"All = "<<endl;
      cout<<"Chi2 = "<<fResolutionAlpha -> GetChisquare()<<endl;
      cout<<"ndof = "<<fResolutionAlpha -> GetNDF()<<endl;
  
      cout<<"qprime full response       = "<<fResolutionAlpha -> GetParameter(1)<<endl;
      cout<<"qprime Error full response = "<<fResolutionAlpha -> GetParError(1)<<endl;
      cout<<"cprime full response       = "<<fResolutionAlpha -> GetParameter(0)<<endl;
      cout<<"cprime Error full response = "<<fResolutionAlpha -> GetParError(0)<<endl;
      cout<<"mprime full response       = "<<fResolutionAlpha -> GetParameter(2)<<endl;
      cout<<"mprime Error full response = "<<fResolutionAlpha -> GetParError(2)<<endl<<endl<<endl;      
    }


    delete fScaleAlpha;
    delete fResolutionAlpha; 
  }
  
};

CScaleResAlpha::~CScaleResAlpha(){
    delete gJetScaleAlpha;
    delete gJetResolutionAlpha;  
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CScaleRes::CScaleRes(){
  
  scalenumber = 0;
  scalenumberError = 0;
  maximum = 0.;
  minimum = 0.;
  for(int i=0; i<nPtBins; i++){
    
    pT[i] = 0;
    pTError[i] = 0;
    scale[i] = 0;
    scaleError[i] = 0;
    resolution[i] = 0;
    resolutionError[i] = 0;
  }
};

void CScaleRes::calculate(int length){

  
  gResolution = new TGraphErrors(length, pT , resolution, pTError, resolutionError);         
  gScale      = new TGraphErrors(length, pT , scale, pTError, scaleError);
  gq          = new TGraphErrors(length, pT , q, pTError, qError);
  gqprime     = new TGraphErrors(length, pT , qprime, pTError, qprimeError);
  
  //====== Scale ==========================================================================================
  fScale = new TF1("fScale","[0]",0,600);  
  gScale -> Fit("fScale","QR");
  scalenumber      = fScale -> GetParameter(0); 
  scalenumberError = fScale -> GetParError(0);
  

  //====== Resolution =====================================================================================
  cout<<endl<<"---------------------------------------------------------------------------------------------------"<<endl;
  fResolution = new TF1("fResolution","TMath::Sqrt(TMath::Sign(1.,[0])*TMath::Power([0]/x,2)+TMath::Power([1],2)*TMath::Power(x,[3]-1)+TMath::Power([2],2))", pT[0]-pTError[0], 600);

  fResolution -> SetParName(0,"N");
  fResolution -> SetParName(1,"S");
  fResolution -> SetParName(2,"C");
  fResolution -> SetParName(3,"m");

  //Set some constraints on the formula (which would lead to unphysical results) and give starting values 
  fResolution -> SetParLimits(fResolution->GetParNumber("m"),-10000.,0.99999999);
  fResolution -> SetParLimits(fResolution->GetParNumber("C"),0.,100000.);
   
  //From Matthias Thesis   
  fResolution->SetParameter("N", -1.45); 
  fResolution->SetParameter("S", 1.034);   
  fResolution->SetParameter("m", -0.085); 
  fResolution->SetParameter("C", 0.05); 
  
  gResolution->Fit("fResolution","QRB");

  maximum = fResolution->GetMaximumX();
  minimum = fResolution->GetMinimumX();

  if(pT[0] < maximum || minimum < 600.){

    cout<<endl<<"Maximum is larger than pT[0]"<<endl;
    cout<<"Maximum is = "<<maximum<<endl;
    cout<<"Minimum is = "<<minimum<<endl;
    cout<<"pT[0] = "<<pT[0]<<endl;
    cout<<"pT[0] - pTError[0] = "<<pT[0]-pTError[0]<<endl;
    
    fResolution -> FixParameter(fResolution->GetParNumber("N"),0.);
    fResolution -> SetParameter("S", 1.034);   
    fResolution -> SetParameter("m", -0.085); 
    fResolution -> SetParameter("C", 0.05); 

    cout<<"Parameter N is set to zero!"<<endl;

    gResolution->Fit("fResolution","QRB");

    maximum = fResolution->GetMaximumX();
    minimum = fResolution->GetMinimumX();
    
    cout<<"Maximum is = "<<maximum<<endl;
    cout<<"Minimum is = "<<minimum<<endl;
    cout<<"pT[0] = "<<pT[0]<<endl<<endl;
  }
   
  N    = fResolution->GetParameter("N");
  NErr = fResolution->GetParError(fResolution->GetParNumber("N"));
  C    = fResolution->GetParameter("C");
  CErr = fResolution->GetParError(fResolution->GetParNumber("C"));
  m    = fResolution->GetParameter("m");  
  mErr = fResolution->GetParError(fResolution->GetParNumber("m"));
  S    = fResolution->GetParameter("S");
  SErr = fResolution->GetParError(fResolution->GetParNumber("S"));

  cout<<"S    = "<<S<<endl;
  cout<<"SErr = "<<SErr<<endl;
  cout<<"m    = "<<m<<endl;
  cout<<"mErr = "<<mErr<<endl;
  cout<<"N    = "<<N<<endl;
  cout<<"NErr = "<<NErr<<endl;  
  cout<<"C    = "<<C<<endl;
  cout<<"CErr = "<<CErr<<endl;
  
  cout<<"---------------------------------------------------------------------------------------------------"<<endl<<endl;
 
};
 

CScaleRes::~CScaleRes(){
  delete gScale;
  delete gResolution;
  delete gq;
  delete gqprime;
  delete fResolution;
  delete fScale;   
};



