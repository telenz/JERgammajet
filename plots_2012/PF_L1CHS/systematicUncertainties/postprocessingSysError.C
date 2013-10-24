// ------------------------------------------------------------------------------------------------
// -------  Script to evaluate all systematic Uncertainties of Jet Energy Resolutions (T.L.)  ------- 
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "../../../CODE/myDeclarations.h"
#include "/afs/naf.desy.de/user/t/telenz/comparison/tdrstyle_mod.C"

//#include "utils.h"
//#include "HistOps.h"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);

int postprocessingSysError(){

  cout<<endl<<endl<<endl<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Postproccess all systematic uncertainties! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl;
  gErrorIgnoreLevel = 1001;

  setTDRStyle(false);

  const TString method  = "RMS99";
  const TString type    = "PFCHS";

  const int nEta =4;
  double eta_bins[5] = {0., 0.5, 1.1, 1.7, 2.3};
  // For looking at different systematic uncertainties independently
  const bool QCD    = true;
  const bool JEC    = true;
  const bool flavor = true;
  const bool PU     = true;
  const bool MC     = true;
  

  TString etaString, filename;   

  TString rootFiles, AuxString;  
  TString JetType = "PFCHS";
  TString Method  = "RMS99";   

  double *ratioEtaBinnedX  = new double[nEta];
  double *ratioEtaBinnedY  = new double[nEta];
  double *ratioEtaBinnedEX = new double[nEta];
  double *ratioEtaBinnedEY = new double[nEta];

  double *ratioEtaBinnedQCDUpY   = new double[nEta];
  double *ratioEtaBinnedQCDDownY = new double[nEta];

  TF1 *QCDuncertainty;

  if(QCD){
    rootFiles   = (TString) "scripts/plotsQCD/FinalErrorsQCD_" + type + (TString) "_" + method + (TString) ".root";
    TFile *_file = TFile::Open(rootFiles);    
    _file->GetObject("function",QCDuncertainty);
  }
  
  for(int eta = 0; eta < nEta; eta++){
    
    //cout<< endl<<endl<<endl<<eta+1<<". eta Bin!!"<<endl;

    // Read the MC and data results 
    rootFiles = (TString) "root_files_FINAL_data/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + JetType + (TString) "_data_" + Method + (TString) ".root";
    TGraphErrors* JERData = readTGraphErrors(rootFiles,"Graph;1","Graph;1");
    rootFiles = (TString) "root_files_FINAL_mc/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + JetType + (TString) "_mc_" + Method + (TString) ".root";
    TGraphErrors* JERMC = readTGraphErrors(rootFiles,"Graph","Graph");
    
    if(eta+1 == 1) etaString = Form("JER for |#eta| < %4.1f",etaBins[eta+1]);
    else           etaString = Form("JER for %4.1f <|#eta|< %4.1f",etaBins[eta+1],etaBins[eta+2]);
 
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // 1.) Calculate the ratio w/o systematic Uncertainties  


    int nData    = JERData->GetN();

    double *dataX  = JERData->GetX();
    double *dataY  = JERData->GetY();
    double *dataEX = JERData->GetEX();
    double *dataEY = JERData->GetEY();
    
    double *mcX  = new double[nData];
    double *mcY  = new double[nData];
    double *mcEX = new double[nData];
    double *mcEY = new double[nData];
    
    double *ratioX  = new double[nData];
    double *ratioY  = new double[nData];
    double *ratioEX = new double[nData];
    double *ratioEY = new double[nData];


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Initialize some stuff for QCD uncertainty
    double *ratioQCDUpY    = new double[nData];
    double *ratioQCDDownY  = new double[nData];    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    int idx = 0;
    for(int i=0; i<nData; i++){

      JERMC   -> GetPoint(idx,mcX[i],mcY[i]);
      mcEX[i] = JERMC -> GetErrorX(idx);
      mcEY[i] = JERMC -> GetErrorY(idx);

      idx += 1;

      if(TMath::Abs(dataX[i]/mcX[i] - 1.) > 0.1){
	i -= 1;
	continue;
      }
     
      ratioX[i]  = 1./2.*(dataX[i] + mcX[i]);
      ratioY[i]  = dataY[i]/mcY[i];
      ratioEX[i] = 1./2.*TMath::Sqrt(TMath::Power(dataEX[i],2)+TMath::Power(mcEX[i],2));
      ratioEY[i] = TMath::Sqrt(TMath::Power((1./mcY[i]),2)*TMath::Power(dataEY[i],2)+TMath::Power((dataY[i]/(TMath::Power(mcY[i],2))),2)*TMath::Power(mcEY[i],2));
      if(QCD){
	// For QCD
	ratioQCDUpY[i]   = ratioY[i]*(1. + QCDuncertainty->Eval(ratioX[i]));
	ratioQCDDownY[i] = ratioY[i]*(1. - QCDuncertainty->Eval(ratioX[i]));
	//cout<<"QCDuncertainty->Eval(ratioX[i]) = "<<QCDuncertainty->Eval(ratioX[i])<<endl;
      }
    }
    
    TGraphErrors *Ratio = new TGraphErrors(nData,ratioX,ratioY,ratioEX,ratioEY);


    // For QCD
    TGraphErrors *QCDUp   = new TGraphErrors(nData,ratioX,ratioQCDUpY,ratioEX,ratioEY);
    TGraphErrors *QCDDown = new TGraphErrors(nData,ratioX,ratioQCDDownY,ratioEX,ratioEY);
    

    if(eta+1 == 1 ) AuxString = Form("Ratio between Data and MC for |#eta| < %4.1f",etaBins[eta+1]);
    else            AuxString = Form("Ratio between Data and MC for %4.1f <|#eta|<%4.1f",etaBins[eta+1],etaBins[eta+2]);
 
    Ratio -> SetTitle(AuxString); 
    Ratio -> GetXaxis() -> SetTitle("Photon pT");
    Ratio -> GetXaxis() -> SetTitleOffset(1.1); 
    Ratio -> GetYaxis() -> SetTitle("Ratio of JER (DATA/MC)");
    Ratio -> GetYaxis() -> SetTitleOffset(1.2);   
    Ratio -> GetXaxis() -> SetLimits(0,600);
    TF1* f1 = new TF1("name","pol0",0,600);   
    Ratio -> Fit("name","QR");
    

    TF1* fitQCDUp  = new TF1("fitQCDUp","pol0",0,600); 
    TF1* fitQCDDown = new TF1("fitQCDDown","pol0",0,600); 
    if(QCD){
      // For QCD
      QCDUp   -> Fit("fitQCDUp","QR");
      QCDDown -> Fit("fitQCDDown","QR");
    }
    
    TLegend *legend  = 0;
    legend = new TLegend(0.55,0.8,0.9,0.9);
    legend -> SetFillColor(0);

    legend -> SetHeader(Form(" %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0)));
    TCanvas *c11 = new TCanvas("c11",AuxString,200,10,500,500);
    c11 -> cd();
    Ratio -> SetMinimum(0.5);
    Ratio -> SetMaximum(2.0);
  
    Ratio  -> Draw("AP"); 
    legend -> Draw("same");
  
    TLatex*  info   = new TLatex();
    info->SetTextFont(132);
    info-> SetNDC();
    info->SetTextSize(0.041); 
    info->DrawLatex(0.22,0.84,Form("#splitline{#chi^{2} = %4.2f}{dof = %i}",f1 -> GetChisquare(),f1 -> GetNDF()));
  
    filename = (TString) "plots/Ratio_Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + type + (TString) "_data_comparison_" + method + (TString) ".pdf";
    c11 -> SaveAs(filename);
    delete c11;
    
    ratioEtaBinnedX[eta]  = (eta_bins[eta+1] + eta_bins[eta])/2.; 
    ratioEtaBinnedY[eta]  = f1 -> GetParameter(0);
    ratioEtaBinnedEX[eta] = 0;
    ratioEtaBinnedEY[eta] = f1->GetParError(0);

    if(QCD){
      ratioEtaBinnedQCDUpY[eta]  = fitQCDUp   -> GetParameter(0);
      ratioEtaBinnedQCDDownY[eta]= fitQCDDown -> GetParameter(0);



      // Some additional stuff for QCD uncertainty
      TCanvas *plotsQCD = new TCanvas("plotsQCD","plotsQCD",200,10,500,500);
      plotsQCD -> cd();

      Ratio -> SetMarkerColor(1);
      Ratio -> SetLineColor(1);
      Ratio -> SetMarkerStyle(20);
      Ratio -> GetFunction("name")->SetLineColor(1);
      QCDUp -> SetMarkerColor(3);
      QCDDown -> SetMarkerColor(3);
      QCDUp  -> SetLineColor(3);
      QCDDown  -> SetLineColor(3);
      QCDUp -> SetMarkerStyle(20);
      QCDDown -> SetMarkerStyle(20);
      QCDUp -> SetMarkerSize(0.8);
      QCDDown -> SetMarkerSize(0.8);
      QCDUp   -> GetFunction("fitQCDUp")->SetLineColor(3);
      QCDDown -> GetFunction("fitQCDDown")->SetLineColor(3);
      Ratio -> Draw("AP");
      QCDUp -> Draw("sameP");
      QCDDown -> Draw("sameP");

      delete legend;
      legend = new TLegend(0.4,0.8,0.9,0.9);
      legend -> SetFillColor(0);
      legend -> SetTextSize(0.033);
      legend -> AddEntry(Ratio,"Central Value","l");
      legend -> AddEntry(QCDUp,Form("Upward variation: + %4.3f",abs(ratioEtaBinnedQCDUpY[eta]/ratioEtaBinnedY[eta]-1.)),"l");
      legend -> AddEntry(QCDDown,Form("Downward variation: - %4.3f",abs(ratioEtaBinnedQCDDownY[eta]/ratioEtaBinnedY[eta]-1.)),"l");
   
      legend -> Draw("same");
      filename = (TString) "plots/plotsQCD_for_" + (long) (eta+1) + (TString) "_bin_"  + type + (TString) "_" + method + (TString) ".pdf";
      plotsQCD -> SaveAs(filename);
      delete plotsQCD;
    }


  }

  TGraphErrors* ratioEtaBinned = new TGraphErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedY,ratioEtaBinnedEX,ratioEtaBinnedEY);
  filename = (TString) "plots/RatioEtaBinned_" + type + (TString) "_" + method + (TString) ".root";
  TFile *f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinned,"Graph");
  f->Close();
  delete f;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 1.) Calculate sys Error from QCD contamination
  //cout<<endl;
    
  double deltaRatioUpQCD[nEta]      = {0.};
  double deltaRatioDownQCD[nEta]    = {0.};
 
  if(QCD){
    
    for(int eta = 0; eta<nEta; eta++){

      deltaRatioUpQCD[eta]     = abs(ratioEtaBinnedQCDUpY[eta]/ratioEtaBinnedY[eta]-1.); 
      deltaRatioDownQCD[eta]   = abs(ratioEtaBinnedQCDDownY[eta]/ratioEtaBinnedY[eta]-1.); 
      
      //cout<<"ratioEtaBinnedQCDDownY[eta]"<<ratioEtaBinnedQCDDownY[eta]<<endl;
      //cout<<"ratioEtaBinnedY[eta]"<<ratioEtaBinnedY[eta]<<endl;
      //cout<<"deltaRatioUpQCD["<<eta<<"] = "<<deltaRatioUpQCD[eta]<<endl;
      //cout<<"deltaRatioDownQCD["<<eta<<"] = "<<deltaRatioDownQCD[eta]<<endl;
    }

  } 
  
  

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 2.) Calculate sys Error from JEC uncertainty (percentage change of MC result)
  //cout<<endl;

  double deltaRatioUpJEC[nEta]      = {0.};
  double deltaRatioDownJEC[nEta]    = {0.};
    
  if(JEC){

    rootFiles                          = (TString) "scripts/plotsJEC/FinalEtaBinnedErrorsJECUp_" + type + (TString) "_" + method + (TString) ".root"; 
    TGraphErrors *JECuncertaintyUp  = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelJECUp       = JECuncertaintyUp -> GetY();

    rootFiles                          = (TString) "scripts/plotsJEC/FinalEtaBinnedErrorsJECLow_" + type + (TString) "_" + method + (TString) ".root"; 
    TGraphErrors *JECuncertaintyLow = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelJECLow      = JECuncertaintyLow -> GetY();
      
    for(int eta = 0; eta<nEta; eta++){

      deltaRatioUpJEC[eta]   = sysRelJECUp[eta];
      deltaRatioDownJEC[eta] = sysRelJECLow[eta];

      //cout<<"deltaRatioUpJEC["<<eta<<"] = "<<deltaRatioUpJEC[eta]<<endl;
      //cout<<"deltaRatioDownJEC["<<eta<<"] = "<<deltaRatioDownJEC[eta]<<endl;     

    }

  }
 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 3.) Calculate sys Error from Flavor uncertainty (percentage change of MC result)
  //cout<<endl;
  
  // Multiply on mc (as symmetric Error)
  // ratioUp  = 1/(1 - delta) * ratio
  // ratioLow = 1/(1 + delta) * ratio

  double deltaRatioUpFlavor[nEta]      = {0.};
  double deltaRatioDownFlavor[nEta]    = {0.};

  if(flavor){

    rootFiles                          = (TString) "scripts/plotsFlavor/FinalEtaBinnedErrorsFlavorUp_" + type + (TString) "_" + method + (TString) ".root"; 
    TGraphErrors *FlavoruncertaintyUp  = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelFlavorUp       = FlavoruncertaintyUp -> GetY();
    
    rootFiles                          = (TString) "scripts/plotsFlavor/FinalEtaBinnedErrorsFlavorLow_" + type + (TString) "_" + method + (TString) ".root"; 
    TGraphErrors *FlavoruncertaintyLow = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelFlavorLow      = FlavoruncertaintyLow -> GetY();
  
    
    for(int eta = 0; eta<nEta; eta++){

      deltaRatioUpFlavor[eta]   = sysRelFlavorUp[eta];
      deltaRatioDownFlavor[eta] = sysRelFlavorLow[eta];

      //cout<<"deltaRatioUpFlavor["<<eta<<"] = "<<deltaRatioUpFlavor[eta]<<endl;
      //cout<<"deltaRatioDownFlavor["<<eta<<"] = "<<deltaRatioDownFlavor[eta]<<endl;
    }

  }
 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 4.) Calculate sys Error from PU uncertainty (percentage change of MC result)
  //cout<<endl;

  double deltaRatioUpPU[nEta]      = {0.};
  double deltaRatioDownPU[nEta]    = {0.};
  
  if(PU){
    
    rootFiles                          = (TString) "scripts/plotsPU/FinalEtaBinnedErrorsPUUp_" + type + (TString) "_" + method + (TString) ".root"; 
    TGraphErrors *PUuncertaintyUp  = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelPUUp       = PUuncertaintyUp -> GetY();

    rootFiles                          = (TString) "scripts/plotsPU/FinalEtaBinnedErrorsPULow_" + type + (TString) "_" + method + (TString) ".root"; 
    TGraphErrors *PUuncertaintyLow = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelPULow      = PUuncertaintyLow -> GetY();
  
    // Multiply on mc (as symmetric Error)
    // ratioUp = 1/(1 - delta) * ratio
    // ratioUp = 1/(1 + delta) * ratio
        
    for(int eta = 0; eta<nEta; eta++){
    
      deltaRatioUpPU[eta]   = sysRelPUUp[eta];
      deltaRatioDownPU[eta] = sysRelPULow[eta];

      //cout<<"deltaRatioUpPU["<<eta<<"] = "<<deltaRatioUpPU[eta]<<endl;
      //cout<<"deltaRatioDownPU["<<eta<<"] = "<<deltaRatioDownPU[eta]<<endl;
    }

  } 
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // 5.) Calculate sys Error from Out-of Cone showering simulation (percentage change of full ratio result)
  //cout<<endl;
  
  double deltaRatioUpMC[nEta]      = {0.};
  double deltaRatioDownMC[nEta]    = {0.};
  

  if(MC){

    rootFiles                   = (TString) "scripts/plotsMC/FinalErrorsMC_" + type + (TString) "_" + method + (TString) ".root";  
    TGraphErrors *MCuncertainty = readTGraphErrors(rootFiles,"graph","Graph");
    double       *sysRelMC      = MCuncertainty -> GetY();
  
    // Percentage change is only in one direction, to take this into account keep deltaRatioDownMC = 0

    for(int eta = 0; eta<nEta; eta++){

      deltaRatioUpMC[eta]   = sysRelMC[eta];
      deltaRatioDownMC[eta] = sysRelMC[eta];

      //cout<<"deltaRatioUpMC["<<eta<<"] = "<<deltaRatioUpMC[eta]<<endl;
      //cout<<"deltaRatioDownMC["<<eta<<"] = "<<deltaRatioDownMC[eta]<<endl;
    }

  } 
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Take all systematic Uncertainties together and plot
  //cout<<endl;

  double *deltaTotalSysUp   = new double[nEta];
  double *deltaTotalSysDown = new double[nEta];
  double *DeltaTotalSysUp   = new double[nEta];
  double *DeltaTotalSysDown = new double[nEta];
  double *DeltaTotalDown = new double[nEta];
  double *DeltaTotalUp = new double[nEta];
  for(int eta = 0; eta<nEta; eta++){

    // Add all systematic Uncertainties in quadrature (delta is relative Uncertainty)
    deltaTotalSysUp[eta]   = sqrt(TMath::Power(deltaRatioUpJEC[eta],2)   + TMath::Power(deltaRatioUpFlavor[eta],2)   + TMath::Power(deltaRatioUpPU[eta],2)   +                                                               TMath::Power(deltaRatioUpMC[eta],2)    + TMath::Power(deltaRatioUpQCD[eta],2));
    deltaTotalSysDown[eta] = sqrt(TMath::Power(deltaRatioDownJEC[eta],2) + TMath::Power(deltaRatioDownFlavor[eta],2) + TMath::Power(deltaRatioDownPU[eta],2) +                                                               TMath::Power(deltaRatioDownMC[eta],2)  + TMath::Power(deltaRatioDownQCD[eta],2));

    // Calculation of the absolute Uncertainty with Delta = ratio * delta
    DeltaTotalSysUp[eta]   = deltaTotalSysUp[eta] * ratioEtaBinnedY[eta];
    DeltaTotalSysDown[eta] = deltaTotalSysDown[eta] * ratioEtaBinnedY[eta];

    // Calculate Systematic plus staistical Uncertainty
    DeltaTotalUp[eta] = sqrt(pow(DeltaTotalSysUp[eta],2) + pow(ratioEtaBinnedEY[eta],2));
    DeltaTotalDown[eta] = sqrt(pow(DeltaTotalSysDown[eta],2) + pow(ratioEtaBinnedEY[eta],2));

    cout<<endl<<"relative: deltaTotalSysUp["<<eta<<"]   = "<<fixed<<setprecision(3)<<deltaTotalSysUp[eta]<<endl;
    cout<<"relative: deltaTotalSysDown["<<eta<<"] = "<<deltaTotalSysDown[eta]<<endl;

    cout<<endl<<"absolute: DeltaTotalSysUp["<<eta<<"]   = "<<DeltaTotalSysUp[eta]<<endl;
    cout<<"absolute: DeltaTotalSysDown["<<eta<<"] = "<<DeltaTotalSysDown[eta]<<endl;


  }

  double ex[nEta] ={0.};

  TGraphAsymmErrors* ratioEtaBinnedSys = new TGraphAsymmErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedY,ex,ex,DeltaTotalSysDown,DeltaTotalSysUp);

  double *TotalSysUp   = new double[nEta];
  double *TotalSysDown = new double[nEta];
 

  for(int i=0; i<nEta; i++){
    TotalSysUp[i]   = ratioEtaBinnedY[i]+DeltaTotalSysUp[i];
    TotalSysDown[i] = ratioEtaBinnedY[i]-DeltaTotalSysDown[i];
  }

  TGraph* ratioSysBorderUp   = new TGraph(nEta, ratioEtaBinnedX, TotalSysUp);
  TGraph* ratioSysBorderDown = new TGraph(nEta, ratioEtaBinnedX, TotalSysDown);

  TGraph* ratioRelativeErrorsUp   = new TGraph(nEta,ratioEtaBinnedX,deltaTotalSysUp);
  TGraph* ratioRelativeErrorsDown = new TGraph(nEta,ratioEtaBinnedX,deltaTotalSysDown);


  TGraphErrors* ratioEtaBinnedStat = new TGraphErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedY,ratioEtaBinnedEX,ratioEtaBinnedEY);
  TGraphAsymmErrors* ratioEtaBinnedStatPlusSys = new TGraphAsymmErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedY,ex,ex,DeltaTotalDown,DeltaTotalUp);
  
  TCanvas *cFinal = new TCanvas("cFinal","cFinal",200,10,500,500);
  cFinal -> cd();  
  
  ratioEtaBinnedSys -> GetYaxis() -> SetTitle("Data/MC ratio (const fit)");
  ratioEtaBinnedSys -> GetXaxis() -> SetTitle("|#eta|");

  if(PU  && flavor  && JEC  && MC && QCD)    etaString  = "All sys. Uncertainties";
  else if(PU  && !flavor && !JEC && !MC && !QCD)  etaString = "Only PU uncert.";
  else if(!PU && flavor  && !JEC && !MC && !QCD)  etaString = "Only flavor uncert.";
  else if(!PU && !flavor && JEC  && !MC && !QCD)  etaString = "Only JEC uncert.";
  else if(!PU && !flavor && !JEC && MC && !QCD)   etaString = "Only Out-of-Cone sim. uncert.";
  else if(!PU && !flavor && !JEC && MC && !QCD)   etaString = "Only Out-of-Cone sim. uncert.";
  else if(!PU && !flavor && !JEC && !MC && QCD)   etaString = "Only QCD uncert.";
  else if(PU && flavor && JEC && !MC && QCD)      etaString = "All besides MC uncertainty.";
  else etaString = "Strange set of systematic uncertainties.";
  cout<<endl<<etaString<<endl<<endl;
  
  ratioEtaBinnedSys -> SetMarkerStyle(20);
  ratioEtaBinnedSys -> SetMarkerSize(1.4);
  ratioEtaBinnedSys -> SetFillColor(kGray);
  ratioEtaBinnedSys -> SetFillStyle(3001);
  ratioEtaBinnedSys -> SetLineColor(kGray);
  ratioEtaBinnedSys -> SetMinimum(0.8);
  ratioEtaBinnedSys -> SetMaximum(1.5);
  ratioEtaBinnedSys -> GetXaxis() -> SetLimits(0., 2.3);
  ratioEtaBinnedSys -> GetXaxis() -> SetNdivisions(6,6,0, "X");
  ratioEtaBinnedSys -> DrawClone("Ae3p");
  
  //ratioEtaBinnedSys -> SetPointError(0, 0., 0., 0., 0.);
  //ratioEtaBinnedSys -> SetPointError(1, 0., 0., 0., 0.);
  //ratioEtaBinnedSys -> SetPointError(2, 0., 0., 0., 0.);
  //ratioEtaBinnedSys -> SetPointError(3, 0., 0., 0., 0.);
  //ratioEtaBinnedSys -> SetPointError(4, 0., 0., 0., 0.);


  ratioEtaBinnedStat -> SetMarkerStyle(20);
  ratioEtaBinnedStat -> SetMarkerSize(1.4);
  ratioEtaBinnedStat -> SetFillColor(kGray);
  ratioEtaBinnedStat -> SetFillStyle(3001);
  ratioEtaBinnedStat -> SetLineColor(1);
  ratioEtaBinnedStat -> Draw("psame");
  
  cmsPrel();

 
  TLatex *infoFinal   = new TLatex();
  infoFinal -> SetTextFont(132);
  infoFinal -> SetNDC();
  infoFinal -> SetTextSize(0.050);
  infoFinal -> DrawLatex(0.2,0.8,etaString);

  filename = (TString) "plots/FinalErrorPlot_" + type + (TString) "_" + method + (TString) ".pdf";
  cFinal -> Print(filename,"pdf");
  filename = (TString) "plots/FinalErrorPlot_" + type + (TString) "_" + method + (TString) ".pdf";
  cFinal -> SaveAs(filename,"pdf");
  delete cFinal;


  filename = (TString) "plots/FinalRelativeErrorsUp_" + type + (TString) "_" + method + (TString) ".root"; 
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioRelativeErrorsUp,"graph");
  f->Close();
  delete f;
  filename = (TString) "plots/FinalRelativeErrorsLow_" + type + (TString) "_" + method + (TString) ".root"; 
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioRelativeErrorsDown,"graph");
  f->Close();
  delete f;

  
  ofstream RelativeErrors;
  RelativeErrors.open("plots/Errors.txt");
  
  RelativeErrors<<"Relative Errors: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<fixed<<setprecision(3)<<(deltaTotalSysDown[i]*100)<<"% / +"<<(deltaTotalSysUp[i]*100)<<"%"<<endl;
  }

  RelativeErrors<<endl<<"Absolute Errors: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(DeltaTotalSysDown[i])<<" / +"<<(DeltaTotalSysUp[i])<<endl;
  }

  RelativeErrors<<endl<<endl<<"Relative Errors JEC: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(deltaRatioDownJEC[i]*100)<<"% / +"<<(deltaRatioUpJEC[i]*100)<<"%"<<endl;
  }
  RelativeErrors<<endl<<endl<<"Relative Errors Flavor: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(deltaRatioDownFlavor[i]*100)<<"% / +"<<(deltaRatioUpFlavor[i]*100)<<"%"<<endl;
  }
  RelativeErrors<<endl<<endl<<"Relative Errors Out-of-Cone showering simulation: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(deltaRatioDownMC[i]*100)<<"% / +"<<(deltaRatioUpMC[i]*100)<<"%"<<endl;
  }
  RelativeErrors<<endl<<endl<<"Relative Errors QCD: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(deltaRatioDownQCD[i]*100)<<"% / +"<<(deltaRatioUpQCD[i]*100)<<"%"<<endl;
  }
  RelativeErrors<<endl<<endl<<"Relative Errors PU reweighing: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(deltaRatioDownPU[i]*100)<<"% / +"<<(deltaRatioUpPU[i]*100)<<"%"<<endl;
  }

  RelativeErrors<<endl<<endl<<"Central values and statistical Uncertainty: "<<endl;
  for(int i=0; i<nEta; i++){
    RelativeErrors<<i+1<<". Eta bin:    "<<"-"<<(ratioEtaBinnedY[i])<<" +/- "<<ratioEtaBinnedEY[i]<<endl;
  }

  RelativeErrors.close();

 
  // Write directly full latex table with systematic and statistical unceratinty

  ofstream latexTable;
  latexTable.open("plots/latexTable.txt");


  latexTable<<"\\renewcommand{\\arraystretch}{2.0}"<<endl;
  latexTable<<"\\begin{center}"<<endl;
  latexTable<<"\\begin{tabular}{ | c | c   c c| }"<<endl;
  latexTable<<"$|\\eta^{\\text{Jet}}|$ & Ratio &  stat.      & sys.  \\\\\\hline"<<endl;
  for(int z=0;z<nEta;z++){
    latexTable<<"$"<<fixed<<setprecision(1)<<etaBins[z]<<" - "<<etaBins[z+1]<<"$ &"<<fixed<<setprecision(3)<<ratioEtaBinnedY[z]<<" & $\\pm "<<ratioEtaBinnedEY[z]<<"$ & $^{+"<<DeltaTotalSysUp[z]<<"}_{-"<<DeltaTotalSysDown[z]<<"}$ \\\\"<<endl;
  }
  latexTable<<"\\hline"<<endl;
  latexTable<<"\\end{tabular}"<<endl;
  latexTable<<"\\end{center}"<<endl<<endl<<endl<<endl<<endl;

 


  
  latexTable<<"\\begin{center}"<<endl;
  latexTable<<"\\begin{tabular}{ l| c | c | c | c |}"<<endl;
  latexTable<<"\\multicolumn{1}{c}{} & \\multicolumn{4}{c}{$|\\eta^{\\text{Jet}}|$}\\\\\\hline"<<endl<<fixed<<setprecision(1);
  for(int z=0;z<nEta;z++) latexTable<<"& \\textbf{"<<etaBins[z]<<" - "<<etaBins[z+1]<<"}";
  latexTable<<"\\\\\\hline"<<endl;
  latexTable<<"\\multirow{2}{*}{\\textbf{Multijet contamination}}";
  for(int z=0;z<nEta;z++) latexTable<<"& $+"<<deltaRatioUpQCD[z]*100<<" \\% $ ";
  latexTable<<"\\\\"<<endl;
  for(int z=0;z<nEta;z++) latexTable<<"& $-"<<deltaRatioDownQCD[z]*100<<" \\% $ ";
  latexTable<<"\\\\\\hline"<<endl;
  latexTable<<"\\multirow{2}{*}{\\textbf{Flavor uncertainty}}";
  for(int z=0;z<nEta;z++) latexTable<<"& $+"<<deltaRatioUpFlavor[z]*100<<" \\% $ ";
  latexTable<<"\\\\"<<endl;
  for(int z=0;z<nEta;z++) latexTable<<"& $"<<deltaRatioDownFlavor[z]*100<<" \\% $ ";
latexTable<<"\\\\\\hline"<<endl;
  latexTable<<"\\multirow{2}{*}{\\textbf{JEC uncertainty}}";
  for(int z=0;z<nEta;z++) latexTable<<"& $+"<<deltaRatioUpJEC[z]*100<<" \\% $ ";
  latexTable<<"\\\\"<<endl;
  for(int z=0;z<nEta;z++) latexTable<<"& $"<<deltaRatioDownJEC[z]*100<<" \\% $ ";
latexTable<<"\\\\\\hline"<<endl;
  latexTable<<"\\multirow{2}{*}{\\textbf{Out-of-Cone showering simulation}}";
  for(int z=0;z<nEta;z++) latexTable<<"& $+"<<deltaRatioUpMC[z]*100<<" \\% $ ";
  latexTable<<"\\\\"<<endl;
  for(int z=0;z<nEta;z++) latexTable<<"& $-"<<deltaRatioDownMC[z]*100<<" \\% $ ";
  latexTable<<"\\\\\\hline"<<endl;
  latexTable<<"\\multirow{2}{*}{\\textbf{PU uncertainty}}";
  for(int z=0;z<nEta;z++) latexTable<<"& $+"<<deltaRatioUpPU[z]*100<<" \\% $ ";
  latexTable<<"\\\\"<<endl;
  for(int z=0;z<nEta;z++) latexTable<<"& $"<<deltaRatioDownPU[z]*100<<" \\% $ ";
  latexTable<<"\\\\\\hline\\hline"<<endl;
  latexTable<<"\\multirow{2}{*}{\\textbf{Total}}";
  for(int z=0;z<nEta;z++) latexTable<<"& $+"<<deltaTotalSysUp[z]*100<<" \\% $ ";
  latexTable<<"\\\\"<<endl;
  for(int z=0;z<nEta;z++) latexTable<<"& $-"<<deltaTotalSysDown[z]*100<<" \\% $ ";
  latexTable<<"\\\\\\hline"<<endl;

  latexTable<<"\\end{tabular}"<<endl;
  latexTable<<"\\end{center}"<<endl;


 latexTable.close();





  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Comparison to 2011 Data 
  cout<<endl; 
  
  TCanvas *cFinal2 = new TCanvas("cFinal2","cFinal2",200,10,1000,1000);
  cFinal2 -> cd();  


  TH1F *Res_2011Final = new TH1F("Data_MC_ratio_2011","", 4, etaBins);
  Res_2011Final->SetBinContent(1, 1.052);
  Res_2011Final->SetBinContent(2, 1.057);
  Res_2011Final->SetBinContent(3, 1.096);
  Res_2011Final->SetBinContent(4, 1.134);
  TGraphAsymmErrors *Res_2011 = new TGraphAsymmErrors(Res_2011Final);
  Res_2011->SetPointError(0, 0., 0., 0.063, 0.062);
  Res_2011->SetPointError(1, 0., 0., 0.057, 0.056);
  Res_2011->SetPointError(2, 0., 0., 0.065, 0.064);
  Res_2011->SetPointError(3, 0., 0., 0.094, 0.092);

  ratioEtaBinnedStatPlusSys -> GetXaxis() -> SetTitle("|#eta|");
  ratioEtaBinnedStatPlusSys -> GetXaxis() -> SetRangeUser(0., 2.3);
  ratioEtaBinnedStatPlusSys -> GetYaxis() -> SetTitle("Data/MC ratio (const fit)");
  ratioEtaBinnedStatPlusSys -> GetYaxis() -> SetRangeUser(0.8, 1.5);

  ratioEtaBinnedSys -> GetXaxis() -> SetTitle("|#eta|");
  ratioEtaBinnedSys -> GetXaxis() -> SetRangeUser(0., 2.3);
  ratioEtaBinnedSys -> GetYaxis() -> SetTitle("Data/MC ratio (const fit)");
  ratioEtaBinnedSys -> GetYaxis() -> SetRangeUser(0.8, 1.5);


  Res_2011 -> GetXaxis() -> SetTitle("|#eta|");
  Res_2011 -> GetXaxis() -> SetLimits(0., 2.3);
  Res_2011 -> GetXaxis() -> SetNdivisions(505, "X");
  Res_2011 -> GetYaxis() -> SetTitle("Data/MC ratio (const fit)");
  Res_2011 -> GetYaxis() -> SetRangeUser(0.8, 1.5);

  ratioEtaBinnedStatPlusSys -> SetMarkerStyle(20); 
  ratioEtaBinnedStatPlusSys -> SetMarkerSize(2.0);
  ratioEtaBinnedStatPlusSys -> SetLineColor(kPink-8);
  ratioEtaBinnedStatPlusSys -> SetMarkerColor(kPink-8);
  ratioEtaBinnedStatPlusSys -> SetFillColor(kPink-8);

  ratioEtaBinnedSys -> SetMarkerStyle(20); 
  ratioEtaBinnedSys -> SetMarkerSize(2.0);
  ratioEtaBinnedSys -> SetLineColor(1);
  ratioEtaBinnedSys -> SetMarkerColor(kPink-8);
  ratioEtaBinnedSys -> SetFillColor(1);

  gStyle->SetHatchesSpacing(2.);
  gROOT->ForceStyle();
  ratioEtaBinnedStatPlusSys -> SetFillStyle(3144);
  ratioEtaBinnedSys -> SetFillStyle(3144);

  Res_2011->SetMarkerStyle(24);
  Res_2011->SetMarkerSize(2.0);
  Res_2011->SetLineColor(1);
  Res_2011->SetFillColor(kGray);
  Res_2011->SetFillStyle(1001);
  Res_2011->SetLineColor(kGray);

  cFinal2->Update();
  Res_2011->DrawClone("Ae3p");
  ratioEtaBinnedStatPlusSys -> Draw("3epsame");
  Res_2011->DrawClone("pXsame");
  Res_2011->SetMarkerSize(1.9);
  Res_2011->DrawClone("pXsame");
  Res_2011->SetMarkerSize(1.7);
  Res_2011->Draw("pXsame");
  ratioEtaBinnedStatPlusSys -> Draw("pXsame");
  
  ratioSysBorderUp -> SetLineColor(12);
  ratioSysBorderUp -> Draw("Lsame");
  ratioSysBorderDown -> SetLineColor(12);
  ratioSysBorderDown -> Draw("Lsame");
    

  //ratioEtaBinnedStatPlusSys -> Draw("e3psame");

  Res_2011->SetPointError(0, 0., 0., 0., 0.);
  Res_2011->SetPointError(1, 0., 0., 0., 0.);
  Res_2011->SetPointError(2, 0., 0., 0., 0.);
  Res_2011->SetPointError(3, 0., 0., 0., 0.);
  Res_2011->SetPointError(4, 0., 0., 0., 0.);

  //Res_2011->Draw("psame");

  //cmsPrel(19.7);
  cmsPrel();

  double m_lumi = 19.7;

  TString infotext = TString::Format("CMS Preliminary, %3.1f fb^{-1}", m_lumi);
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetX(0.22);
  text1->SetTextFont(42);
  infotext = TString::Format("#sqrt{s} = 8 TeV");
  TLatex *text2 = new TLatex(3.5, 24, infotext);
  text2->SetNDC();
  text2->SetX(0.22);
  text2->SetTextFont(42);

  text1->SetTextSize(0.040);
  text1->SetTextAlign(11);
  text1->SetY(0.96);
  text1->SetX(0.15);

  text2->SetTextSize(0.040);
  text2->SetTextAlign(31);
  text2->SetY(0.96);
  text2->SetX(0.95);

  text1->Draw("same");
  text2->Draw("same");
  

  TLegend *leg = new TLegend(0.20, 0.70, 0.40, 0.90);
  leg->SetBorderSize(0);
  // leg->SetBorderMode(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.040);
  
  leg->AddEntry(Res_2011,"2011 (total error)", "plf");
  leg->AddEntry(ratioEtaBinnedStatPlusSys,"2012 (total error)", "plf");
  leg->AddEntry(ratioSysBorderUp,"systematic error (2012)", "l");
 
   
  leg->Draw("same");
  
  cFinal2->Print("plots/resultsComparisonFINAL.pdf","pdf");


  return 0;

}













//! Read TGraphErrors from file
// -------------------------------------------------------------------------------------
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraphErrors *g = 0;
  file.GetObject(gName,g);
  if( g ) {
    if( newGName.Length() ) g->SetName(newGName);
  } else {
    std::cerr << "ERROR in FileOps::readTGraph: TGraphAsymmErrors with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    file.Close();
    exit(-1);
  }
  file.Close();
    
  return g;
} 

//! Read TGraph from file
// -------------------------------------------------------------------------------------
TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraph *g = 0;
  file.GetObject(gName,g);
  if( g ) {
    if( newGName.Length() ) {
    }  
    g->SetName(newGName);
  } else {
    std::cerr << "ERROR in FileOps::readTGraph: TGraph with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    file.Close();
    exit(-1);
  }
  file.Close();
    
  return g;
}

//! Read TGraphAsymmErrors from file
// -------------------------------------------------------------------------------------
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraphAsymmErrors *g = 0;
  file.GetObject(gName,g);
  if( g ) {
    if( newGName.Length() ) g->SetName(newGName);
  } else {
    std::cerr << "ERROR in FileOps::readTGraph: TGraphAsymmErrors with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
    file.Close();
    exit(-1);
  }
  file.Close();
    
  return g;
}









