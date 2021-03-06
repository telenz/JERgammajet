// $Id: evalJECUncertainty.C,v 1.2 2013/07/10 08:54:28 telenz Exp $

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "../../../../CODE/myDeclarations.h"
#include "../../../../scripts/plotStyle.h"
#include "../../../../CODE/myFunctions.h"
#include "../../../../CODE/myFunctions.C"


int evalJECUncertainty(){

  cout<<endl<<endl<<endl<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for JEC uncertainty is executed! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl;
  gErrorIgnoreLevel = 1001;

  TeresaPlottingStyle::init();

  const int nEta = 4;

  const TString method = "RMS99";
  const TString type   = "PFCHS";
 
  TString etaString = "Uncertainty on JEC";
  
  TString tot_filename, AuxString, fitName;;
  
  double finalErrorsUpY[nEta][nPtBins]  = {{0}};
  double finalErrorsLowY[nEta][nPtBins] = {{0}};
  double finalErrorsUpX[nEta][nPtBins]  = {{0}};
  double finalErrorsLowX[nEta][nPtBins] = {{0}};
  int nCount[nEtaBins] = {0};

  TString rootFile[3]; 
  TString pathName[3];
 

  /*  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Calculate Correlation with the help of the control plot hPhoton1Pt_PFCHS_mc.root
  
  //rootFile[0]  = (TString) "woSysUncer_woTriggerWithPixelSeed/hPhoton1Pt_" + type + (TString) "_mc.root";
  //rootFile[1]  = (TString) "JECUncertainties/downwardVariation_woTriggerWithPixelSeed/hPhoton1Pt_" + type + (TString) "_mc.root";
  //rootFile[2]  = (TString) "JECUncertainties/upwardVariation_woTriggerWithPixelSeed/hPhoton1Pt_" + type + (TString) "_mc.root"; 

  rootFile[0]  = (TString) "root_files_woTriggerPUeq1_2ndJet/hPhoton1Pt_" + type + (TString) "_mc.root";
  rootFile[1]  = (TString) "JECUncertainties/downwardVariation_woTriggerPUeq1_2ndJet/hPhoton1Pt_" + type + (TString) "_mc.root";
  rootFile[2]  = (TString) "JECUncertainties/upwardVariation_woTriggerPUeq1_2ndJet/hPhoton1Pt_" + type + (TString) "_mc.root";
  
  TFile* file[3];
  TH1D* histo[3];
  
  for(int i =0; i<3; i++){
  file[i] = TFile::Open(rootFile[i]);
  file[i]->GetObject("histo",histo[i]);    
  }
   
  double correlationLow = 1.;
  double correlationUp  = 1.;
  if(histo[1]->Integral() < histo[0]->Integral()) correlationLow = sqrt(histo[1]->Integral()/histo[0]->Integral());
  else                                            correlationLow = sqrt(histo[0]->Integral()/histo[1]->Integral());
  if(histo[2]->Integral() < histo[0]->Integral()) correlationUp  = sqrt(histo[2]->Integral()/histo[0]->Integral());
  else                                            correlationUp  = sqrt(histo[0]->Integral()/histo[2]->Integral());
  */
  double correlationLow = 1.;
  double correlationUp  = 1.;
  
  cout<<endl<<"Correlation between downward variation and without variation = "<<correlationLow<<endl;
  cout<<      "Correlation between upward variation and without variation   = "<<correlationUp<<endl<<endl<<endl;
  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
  pathName[0]  = (TString) "root_files_JECUncertainty/centralValue/";
  pathName[1]  = (TString) "root_files_JECUncertainty/downwardVariation/";
  pathName[2]  = (TString) "root_files_JECUncertainty/upwardVariation/";

  cout<<"root files from following folders:"<<endl<<pathName[0]<<endl<<pathName[1]<<endl<<pathName[2]<<endl<<endl<<endl;

  for(int eta = 0; eta<nEta; eta++){
 
    for(int i=0; i<3; i++) rootFile[i]  = pathName[i] + (TString) "Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_intrinsic_" + type + (TString) "_mc_" + method + (TString) ".root";
  
    TMultiGraph* mg = new TMultiGraph();
    etaString       = "Uncertainty on JEC";
    mg             -> SetTitle(etaString);
    

    TCanvas *c = new TCanvas("c",etaString,200,10,800,800);
    c -> cd();
  
    double maximum = 0;
   
    TGraphErrors* graph[3]; 
    for(int i =0; i<3; i++){
      graph[i] = readTGraphErrors(rootFile[i],"Graph","Graph");    
      if(graph[i]->GetYaxis()->GetXmax() > maximum) maximum = graph[i]->GetYaxis()->GetXmax();
    }
  
    graph[0] -> SetMarkerColor(8);
    graph[1] -> SetMarkerColor(9);
    graph[2] -> SetMarkerColor(1);
 
    graph[0] -> SetLineColor(8);
    graph[1] -> SetLineColor(9);
    graph[2] -> SetLineColor(1);

    graph[0] -> GetFunction("fResolution") -> SetLineColor(8);
    graph[1] -> GetFunction("fResolution") -> SetLineColor(9);
    graph[2] -> GetFunction("fResolution") -> SetLineColor(1);
  
 
    TLegend *legend  = new TLegend(0.5,0.8,0.9,0.9);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(graph[0],"Central Value","l");
    legend -> AddEntry(graph[1],"- #Delta JEC","l");
    legend -> AddEntry(graph[2],"+ #Delta JEC","l");
  
    for(int i =0; i<3; i++) mg->Add(graph[i]);  
      
    mg->Draw("AP");
    
    mg -> GetYaxis() -> SetTitle("JER");
    mg -> SetMinimum(0.004);
    mg -> SetMaximum(0.12);   
    mg -> GetXaxis() -> SetLimits(0,600);
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]"); 
    legend->Draw("same");

    TLatex*  info   = new TLatex();
    
    info-> SetNDC();
    AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta],etaBins[eta+1]);
    info->DrawLatex(0.6,0.7,AuxString);


    tot_filename = (TString) "plotsJEC/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_JECUncertainty_" + method + (TString) ".pdf";
    c -> SaveAs(tot_filename);
 
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // 2.) Relative uncertainty to MC -> MC*(1 +- Delta)
    double *dataY, *dataLowY, *dataUpY, *dataX, *dataLowX, *dataUpX, *dataEY, *dataLowEY, *dataUpEY, *dataEX;

    int nData = graph[0]->GetN();
  
    if(nData > graph[1]->GetN()) nData = nData - (nData - graph[1]->GetN()); 
    if(nData > graph[2]->GetN()) nData = nData - (nData - graph[2]->GetN());  

    double *errorUpY  = new double[nData];
    double *errorLowY = new double[nData];
    double *errorUpEY = new double[nData];
    double *errorLowEY= new double[nData];

    dataY     = graph[0] -> GetY();
    dataLowY  = graph[1] -> GetY();
    dataUpY   = graph[2] -> GetY();
    dataX     = graph[0] -> GetX();
    dataLowX  = graph[1] -> GetX();
    dataUpX   = graph[2] -> GetX();
    dataEY    = graph[0] -> GetEY();
    dataLowEY = graph[1] -> GetEY();
    dataUpEY  = graph[2] -> GetEY();
    dataEX    = graph[0] -> GetEX();

    int idxLow     = 0;
    int idxUp      = 0;
    int idx        = 0;
    int countNData = 0;

    while(countNData < nData){


      if(TMath::Abs(dataX[idx]/dataLowX[idxLow] - 1.) > 0.05){
	if(dataX[idx] > dataLowX[idxLow]) idxLow += 1;
	else idx += 1;
	continue;
      }
      if(TMath::Abs(dataX[idx]/dataUpX[idxUp] - 1.) > 0.05){
	if(dataX[idx] > dataUpX[idxUp]) idxUp += 1;
	else idx += 1;
	continue;
      }

      if(countNData == 0 || countNData == nData -1){
	//cout<<"dataUpX["<<idxUp<<"] = "<<dataUpX[idxUp]<<endl;
	//cout<<"dataLowX["<<idxLow<<"] = "<<dataLowX[idxLow]<<endl;
	//cout<<"dataX["<<idx<<"] = "<<dataX[idx]<<endl;
      }

      errorUpY[countNData]   = dataUpY[idxUp]/dataY[idx] - 1.;
      errorLowY[countNData]  = dataLowY[idxLow]/dataY[idx] - 1;

      errorUpEY[countNData]  = sqrt(pow(1./dataY[idx],2)*pow(dataUpEY[idxUp],2) + pow(dataUpY[idxUp]/pow(dataY[idx],2),2)*pow(dataEY[idx],2) - correlationUp*2.* dataUpY[idxUp]/pow(dataY[idx],3)*dataUpEY[idxUp]*dataEY[idx]);
      errorLowEY[countNData] = sqrt(pow(1./dataY[idx],2)*pow(dataLowEY[idxLow],2) + pow(dataLowY[idxLow]/pow(dataY[idx],2),2)*pow(dataEY[idx],2) - correlationLow*2.0 * dataLowY[idxLow]/pow(dataY[idx],3)*dataLowEY[idxLow]*dataEY[idx]);

      dataX[countNData]  = dataX[idx];
      dataEX[countNData] = dataEX[idx];

      countNData += 1;
      idx        += 1;
      idxUp      += 1;
      idxLow     += 1;
    }


    TGraphErrors *plotUp  = new TGraphErrors(nData,dataX,errorUpY,dataEX,errorUpEY);
    TGraphErrors *plotLow = new TGraphErrors(nData,dataX,errorLowY,dataEX,errorLowEY);

    //TF1 *fitLineLow  =  new TF1("fitLineLow","expo",0,600);
    //TF1 *fitLineUp   =  new TF1("fitLineUp","expo",0,600);
    //plotUp -> Fit("fitLineUp","QR");
    //plotLow -> Fit("fitLineLow","QR");

    
    double *xUp   = plotUp  -> GetX();
    double *xLow  = plotLow -> GetX();
    double *yUp   = plotUp  -> GetY();
    double *yLow  = plotLow -> GetY();

    int idxPtBins = 0;
    nCount[eta] = 0;

    for(int pt=0;pt<nPtBins;pt++){
      
      if(idxPtBins>nPtBins) break;

      if(!(xUp[pt]<ptBins[idxPtBins+1] && xUp[pt]>ptBins[idxPtBins])){
	idxPtBins += 1;
	pt -= 1;
	continue;
      }

      finalErrorsUpX[eta][pt]  = xUp[pt];
      finalErrorsUpY[eta][pt]  = yUp[pt];
      finalErrorsLowX[eta][pt] = xLow[pt];
      finalErrorsLowY[eta][pt] = yLow[pt];

      nCount[eta] += 1;
    }

    TF1 *Line   =  new TF1("Line","pol0",0.,600.);
    Line->SetParameter(0,0);
  
    TCanvas *c1 = new TCanvas("c1","c1",200,10,800,800);
    c1 -> cd();
    c1 -> SetBottomMargin(0.14);
    c1 -> SetLeftMargin(0.14);

    plotUp -> SetMarkerStyle(20);  
    plotUp -> SetMarkerColor(1);  
    plotUp -> SetLineColor(1); 

    plotLow -> SetMarkerStyle(24);  
    plotLow -> SetMarkerColor(9);  
    plotLow -> SetLineColor(9); 
 
    TLegend *legend1  = new TLegend(0.5,0.8,0.9,0.9);
    legend1 -> SetTextSize(0.033);
    legend1 -> AddEntry(plotLow,"- #Delta JEC","p");
    legend1 -> AddEntry(plotUp,"+ #Delta JEC","p");
 
    mg = new TMultiGraph();
    mg->Add(plotLow);  
    mg->Add(plotUp);
    //-----

    TLatex* info1   = new TLatex();
    info1->SetTextFont(132);
    info1-> SetNDC();
    info1->DrawLatex(0.6,0.7,AuxString);
    
    //---

    etaString = "Relative JEC Uncertainty";
    mg -> SetTitle(etaString);  
  
    mg->Draw("AP");

    mg -> GetYaxis() -> SetTitle("JER_{#pm #Delta JEC} /JER -1");
    mg -> SetMinimum(-0.1);
    mg -> SetMaximum(0.1); 
    mg -> GetXaxis() -> SetLimits(0,600);  
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");

    legend1->Draw("same");

    delete info1;
    info1   = new TLatex();
    info1-> SetNDC();

    info1->DrawLatex(0.6,0.7,AuxString);

    Line->Draw("same");

    tot_filename = (TString) "plotsJEC/Relative_Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_JECUncertainty_" + method + (TString) ".pdf";
  
    c1 -> SaveAs(tot_filename);

    delete mg;  
    delete c1;
    delete c;

  }

  // Save relative uncertainties for every eta bin in another root-file
  double eta[nEta] = {1.};
  double etaError[nEta] = {0.};
  for(int i =0; i<nEta-1; i++) eta[i+1] = eta[i]+1.;


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Move MC points up and down with the errors and fit again the data/MC_up and data/MC_down ratios -> Get errors for every eta Bin

  TString rootFiles;
  TString JetType = "PFCHS";

  double *ratioEtaBinnedX     = new double[nEta];
  double *ratioEtaBinnedY     = new double[nEta];
  double *ratioEtaBinnedUpY   = new double[nEta];
  double *ratioEtaBinnedLowY  = new double[nEta];
  double *ratioEtaBinnedEX    = new double[nEta];
  double *ratioEtaBinnedEY    = new double[nEta];
  double *ratioEtaBinnedUpEY  = new double[nEta];
  double *ratioEtaBinnedLowEY = new double[nEta];

  for(int iEta = 0; iEta < nEta; iEta++){
    
    // Read the MC and data results 
    rootFiles = (TString) "../root_files_FINAL_data/Resolution_for_" + (long) (iEta+1) + (TString) "_eta_bin_" + JetType + (TString) "_data_" + method + (TString) ".root";
    TGraphErrors* JERData = readTGraphErrors(rootFiles,"Graph;1","Graph;1");
    rootFiles = (TString) "../root_files_FINAL_mc/Resolution_for_" + (long) (iEta+1) + (TString) "_eta_bin_" + JetType + (TString) "_mc_" + method + (TString) ".root";
    TGraphErrors* JERMC = readTGraphErrors(rootFiles,"Graph","Graph");   
    if(iEta+1 == 1) etaString = Form("JER for |#eta| < %4.1f",etaBins[iEta+1]);
    else           etaString = Form("JER for %4.1f <|#eta|< %4.1f",etaBins[iEta+1],etaBins[iEta+2]);
 
    int nData    = JERData->GetN();
    
    double *dataX  = JERData->GetX();
    double *dataY  = JERData->GetY();
    double *dataEX = JERData->GetEX();
    double *dataEY = JERData->GetEY();
    
    double *mcX     = new double[nData];
    double *mcY     = new double[nData];
    double *dataUpY   = new double[nData];
    double *dataLowY  = new double[nData];
    double *mcEX    = new double[nData];
    double *mcEY    = new double[nData];
        
    double *ratioX     = new double[nData];
    double *ratioY     = new double[nData];
    double *ratioUpY   = new double[nData];
    double *ratioLowY  = new double[nData];
    double *ratioEX    = new double[nData];
    double *ratioEY    = new double[nData];
     
    int idx     = 0;
    int idxUp   = 0;
    int idxLow  = 0;
    int idxData = 0;
    int idxMC   = 0;
    
    while(idx<nData){

      if(idxMC>nPtBins || idxUp>nPtBins || idxLow>nPtBins){
	//cout<<"idxMC = "<<idxMC<<endl;
	//cout<<"idxLow = "<<idxLow<<endl;
	//cout<<"idxUp = "<<idxUp<<endl;
	//cout<<"nData = "<<nData<<endl;
	//cout<<"nPtBins = "<<nPtBins<<endl;
	cout<<"something wrong here"<<endl;
	break;
      }

      if(abs(dataX[idxData]/finalErrorsUpX[iEta][idxUp]-1.)>0.05){
	idxUp += 1;
	continue;
      }
      else{
	dataUpY[idxData] = dataY[idxData]*(1.+finalErrorsUpY[iEta][idxUp]);
      }

      if(abs(dataX[idxData]/finalErrorsLowX[iEta][idxLow]-1.)>0.05){
	idxLow += 1;
	continue;
      }
      else{
	dataLowY[idxData] = dataY[idxData]*(1.+finalErrorsLowY[iEta][idxLow]);
      }


      JERMC       -> GetPoint(idxMC,mcX[0],mcY[0]);
      mcEX[0]   = JERMC -> GetErrorX(idxMC);
      mcEY[0]   = JERMC -> GetErrorY(idxMC);

      if(abs(dataX[idxData]/mcX[0] - 1.) > 0.05){

	if(dataX[idxData]<mcX[0]){
	  idxData += 1;
	  nData -= 1;
	  continue; 
	}
	else{
	  idxMC  += 1;
	  //idxUp  += 1;
	  //idxLow += 1;
	  continue;
	}
      }
          
      ratioX[idx]  = 1./2.*(dataX[idxData] + mcX[0]);
      ratioY[idx]  = dataY[idxData]/mcY[0];
      ratioUpY[idx]   = dataLowY[idxData]/mcY[0];
      ratioLowY[idx]  = dataUpY[idxData]/mcY[0];
      ratioEX[idx] = 1./2.*TMath::Sqrt(TMath::Power(dataEX[idxData],2)+TMath::Power(mcEX[0],2));
      ratioEY[idx] = TMath::Sqrt(TMath::Power((1./mcY[0]),2)*TMath::Power(dataEY[idxData],2)+TMath::Power((dataY[idxData]/(TMath::Power(mcY[0],2))),2)*TMath::Power(mcEY[0],2));


      idxData += 1;
      idxMC   += 1;
      idxUp   += 1;
      idxLow  += 1;
      idx     += 1;
    }
    
    TGraphErrors *Ratio    = new TGraphErrors(nData,ratioX,ratioY,ratioEX,ratioEY);
    TGraphErrors *RatioUp  = new TGraphErrors(nData,ratioX,ratioUpY,ratioEX,ratioEY);
    TGraphErrors *RatioLow = new TGraphErrors(nData,ratioX,ratioLowY,ratioEX,ratioEY);

    if(iEta+1 == 1 ) AuxString = Form("Ratio between Data and MC for |#eta| < %4.1f",etaBins[iEta+1]);
    else             AuxString = Form("Ratio between Data and MC for %4.1f <|#eta|<%4.1f",etaBins[iEta+1],etaBins[iEta+2]);
 
    Ratio -> SetTitle(AuxString); 
    Ratio -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
    Ratio -> GetXaxis() -> SetTitleOffset(1.1); 
    Ratio -> GetYaxis() -> SetTitle("Ratio of JER (DATA/MC)");
    Ratio -> GetYaxis() -> SetTitleOffset(1.2);   
    Ratio -> GetXaxis() -> SetLimits(0,600);
    TF1* f1 = new TF1("name","pol0",0,600);  
    f1->SetLineColor(2); 
    Ratio->SetMarkerColor(2); 
    Ratio->SetLineColor(2); 
    Ratio -> Fit("name","QR");
    TF1* f1Up = new TF1("nameUp","pol0",0,600);   
    f1Up->SetLineColor(9); 
    RatioUp->SetMarkerColor(9);
    RatioUp->SetLineColor(9);
    RatioUp -> Fit("nameUp","QR");
    TF1* f1Low = new TF1("nameLow","pol0",0,600);   
    f1Low->SetLineColor(1);
    RatioLow->SetMarkerColor(1); 
    RatioLow->SetLineColor(1); 
    RatioLow -> Fit("nameLow","QR");
    

    TLegend *legend  = 0;
    legend = new TLegend(0.50,0.7,0.9,0.9);
    legend -> SetFillColor(0);

    legend -> AddEntry(Ratio,Form("Ratio =  %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0)),"pl");
    legend -> AddEntry(RatioUp,Form("RatioUp =  %4.3f #pm %4.3f", f1Up -> GetParameter(0), f1Up->GetParError(0)),"pl");
    legend -> AddEntry(RatioLow,Form("RatioDown =  %4.3f #pm %4.3f", f1Low -> GetParameter(0), f1Low->GetParError(0)),"pl");


   
    TCanvas *c11 = new TCanvas("c11",AuxString,200,10,500,500);
    c11 -> cd();
    Ratio -> SetMinimum(0.85);
    Ratio -> SetMaximum(1.25);
  
    Ratio  -> Draw("AP"); 
    RatioUp  -> Draw("Psame");
    RatioLow  -> Draw("Psame");
    legend -> Draw("same");

    TLatex*  info   = new TLatex();
    info-> SetNDC();
    AuxString = Form("#splitline{+%4.1f %%}{%4.1f %%}",(f1Up->GetParameter(0)/f1->GetParameter(0)-1)*100,(f1Low->GetParameter(0)/f1->GetParameter(0)-1)*100);
    info->DrawLatex(0.6,0.30,AuxString);


    
    TString filename = (TString) "plotsJEC/Ratios_upper_lower_Error_for_" + (long) (iEta+1) + (TString) "_eta_bin_" + type + (TString) "_" + method + (TString) ".pdf";
    c11 -> SaveAs(filename);
    delete c11;
 
    
    ratioEtaBinnedX[iEta]     = (etaBins[iEta+1] + etaBins[iEta])/2.; 
    ratioEtaBinnedY[iEta]     = f1 -> GetParameter(0);
    ratioEtaBinnedUpY[iEta]   = f1Up -> GetParameter(0)/f1 -> GetParameter(0)-1.;
    ratioEtaBinnedLowY[iEta]  = f1Low -> GetParameter(0)/f1 -> GetParameter(0)-1.;
    ratioEtaBinnedEX[iEta]    = 0;
    ratioEtaBinnedEY[iEta]    = f1->GetParError(0);
    ratioEtaBinnedUpEY[iEta]  = f1Up->GetParError(0);
    ratioEtaBinnedLowEY[iEta] = f1Low->GetParError(0);

    
    cout<<"upper rel. Error "<<iEta<<". eta Bin = "<<ratioEtaBinnedUpY[iEta]<<endl;
    cout<<"lower rel. Error "<<iEta<<". eta Bin = "<<ratioEtaBinnedLowY[iEta]<<endl<<endl;
  }

  TGraphErrors* ratioEtaBinned    = new TGraphErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedY,ratioEtaBinnedEX,ratioEtaBinnedEY);
  TGraphErrors* ratioEtaBinnedUp  = new TGraphErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedUpY,ratioEtaBinnedEX,ratioEtaBinnedUpEY);
  TGraphErrors* ratioEtaBinnedLow = new TGraphErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedLowY,ratioEtaBinnedEX,ratioEtaBinnedLowEY);

  TString filename = (TString) "plotsJEC/RatioEtaBinned_" + type + (TString) "_" + method + (TString) ".root";
  TFile* f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinned,"Graph");
  f->Close();
  delete f;

  filename = (TString) "plotsJEC/FinalEtaBinnedErrorsJECUp_" + type + (TString) "_"  + method + (TString) ".root"; 
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinnedUp,"graph");
  f->Close();
  delete f;

  filename = (TString) "plotsJEC/FinalEtaBinnedErrorsJECLow_" + type + (TString) "_"  + method + (TString) ".root"; 
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinnedLow,"graph");
  f->Close();
  delete f;

  return 0;
}






