// PU Uncertainty evaluation

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
#include "../../../../CODE/myDeclarations.h"
#include "../../../../scripts/plotStyle.h"
#include "../../../../CODE/myFunctions.h"
#include "../../../../CODE/myFunctions.C"


int evalPUUncertainty(){

  cout<<endl<<endl<<endl<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for PU uncertainty is executed! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl;
  gErrorIgnoreLevel = 1001;

  TeresaPlottingStyle::init();

  const int nEta = 4;

  const TString method = "RMS99";
  const TString type   = "PFCHS";
 
  TString etaString = "PU Uncertainty";
  
  TString tot_filename, AuxString, fitName;;
  
  double finalErrors[nEta]  = {0};
  double finalErrorsE[nEta] = {0};

  TString rootFile[3];
  TString pathName[3];
 
  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  double correlationLow = 1.;
  double correlationUp  = 1.;
  
  cout<<endl<<"Correlation between downward variation and without variation = "<<correlationLow<<endl;
  cout<<      "Correlation between upward variation and without variation   = "<<correlationUp<<endl<<endl<<endl;
  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
  pathName[0]  = (TString) "root_files_PUUncertainty/centralValue/";
  pathName[1]  = (TString) "root_files_PUUncertainty/downwardVariation/";
  pathName[2]  = (TString) "root_files_PUUncertainty/upwardVariation/";

  cout<<"root files from following folders:"<<endl<<pathName[0]<<endl<<pathName[1]<<endl<<pathName[2]<<endl<<endl<<endl;


  for(int eta = 1; eta<nEta+1; eta++){

    for(int i=0; i<3; i++) rootFile[i]  = pathName[i] + (TString) "Resolution_for_" + (long) eta + (TString) "_eta_bin_" + type + (TString) "_mc_" + method + (TString) ".root";
   
    TMultiGraph* mg = new TMultiGraph();
    etaString       = "PU Uncertainty";
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
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.042);
    legend -> AddEntry(graph[0],"MBX = 69.4 mb","l");
    legend -> AddEntry(graph[1],"MBX = 65.8 mb","l");
    legend -> AddEntry(graph[2],"MBX = 73.0 mb","l");
  
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
    AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta-1],etaBins[eta]);
    info->DrawLatex(0.6,0.7,AuxString);


    tot_filename = (TString) "plotsPU/Resolution_for_" + (long) eta + (TString) "_eta_bin_PUUncertainty_" + method + (TString) ".pdf";
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

    // Search for the interval where 68% of all point are included
    double* arraySort = new double[2*nData];
    float auxSort     = 10000000.;

  
    for(int i=0; i<nData; i++){
      arraySort[i]       = abs(errorUpY[i]);
      arraySort[i+nData] = abs(errorLowY[i]);
    }

    for(int j=0; j<2*nData; j++){
      for(int i=0; i<2*nData-1; i++){
	if(arraySort[i] > arraySort[i+1]){
	  auxSort        = arraySort[i];
	  arraySort[i]   = arraySort[i+1];
	  arraySort[i+1] = auxSort;
	}
      }
    }

    TF1* sigma1Interval = new TF1("sigma1Interval","pol0",0,600);
    sigma1Interval->SetLineColor(8);

    // 0.5 are added to have the right rounding to an integer
    int interval68 = (int) (2.*nData*0.6827 + 0.5);
    int interval95 = (int) (2.*nData*0.9545 + 0.5);
    //cout<<endl<<"2.*nData*0.6827 = "<<2.*nData*0.6827<<endl;
    //cout<<"2.*nData*0.9545 = "<<2.*nData*0.9545<<endl;
    //cout<<"interval95 = "<<interval95<<endl;
    //cout<<"interval68 = "<<interval68<<endl<<endl<<endl;
  
  

    if(2*arraySort[interval68-1]<arraySort[interval95-1]){
      sigma1Interval->SetParameter(0,arraySort[interval95-1]/2.);
      //cout<<"arraySort["<<interval68-1<<"]    = "<<arraySort[interval68-1]<<endl;
      //cout<<"arraySort["<<interval95-1<<"]/2. = "<<arraySort[interval95-1]/2.<<endl;
      //cout<<"arraySort["<<2*nData-1<<"] = "<<arraySort[2*nData -1]<<endl;
      //cout<<"arraySort["<<2*nData-2<<"] = "<<arraySort[2*nData -2]<<endl;
    }
    else sigma1Interval->SetParameter(0,arraySort[interval68-1]);
    sigma1Interval->SetParameter(0,arraySort[interval68-1]);
    cout<<"sigma1Interval->GetParameter(0) = "<<sigma1Interval->GetParameter(0)<<endl<<endl;
    TF1* sigma1Intervaldown = new TF1("sigma1Intervaldown","pol0",0,600);
    sigma1Intervaldown->SetParameter(0,-sigma1Interval->GetParameter(0));
    sigma1Intervaldown->SetLineColor(8);
  
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
    legend1 -> SetFillColor(0);
    legend1 -> SetTextSize(0.042);
    legend1 -> AddEntry(plotLow,"MBX = 65.8 mb","p");
    legend1 -> AddEntry(plotUp,"MBX = 73.0 mb","p");

 
    mg = new TMultiGraph();
    mg->Add(plotLow);  
    mg->Add(plotUp);
  

    etaString = "Relative PU Uncertainty";
    mg -> SetTitle(etaString);  
  
    mg->Draw("AP");
    sigma1Interval->Draw("same");
    sigma1Intervaldown->Draw("same");    

    mg -> GetYaxis() -> SetTitle("JER_{MBX = 73.0/65.8} /JER_{MBX = 69.4}");
    mg -> SetMinimum(-0.05);
    mg -> SetMaximum(0.05); 
    mg -> GetXaxis() -> SetLimits(0,600);  
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  
    legend1->Draw("same");

    TLatex*  info1   = new TLatex();
    info1-> SetNDC();

    info1->DrawLatex(0.6,0.7,AuxString);
    TString text = "#scale[0.9]{68% of points are included in Interval:} " ;
    info1->DrawLatex(0.2,0.3,text);
    text = (TString) Form("#scale[0.9]{#delta^{PU} = #pm %4.3f}",sigma1Interval->GetParameter(0)) ;
    info1->DrawLatex(0.7,0.25,text);

    finalErrors[eta-1] = sigma1Interval->GetParameter(0);
    finalErrorsE[eta-1] = sigma1Interval->GetParError(0);



    tot_filename = (TString) "plotsPU/Relative_Resolution_for_" + (long) eta + (TString) "_eta_bin_PUUncertainty_" + method + (TString) ".pdf";
  
    c1 -> SaveAs(tot_filename);

    delete mg;  
    delete c1;
    delete c;

  }

  // Save relative uncertainties for every eta bin in another root-file
  double eta[nEta] = {1.};
  double etaError[nEta] = {0.};
  for(int i =0; i<nEta-1; i++) eta[i+1] = eta[i]+1.;

  
  TGraphErrors* finalErrorsPU = new TGraphErrors(nEta,eta,finalErrors,etaError,finalErrorsE);
  finalErrorsPU -> SetMarkerStyle(20);
  finalErrorsPU -> SetTitle("Final relative Erros (PU)");
  finalErrorsPU -> GetXaxis() -> SetTitle("#eta^{Jet}");
  finalErrorsPU -> GetYaxis() -> SetTitle("JER_{MBX = 73.0/65.8} /JER_{MBX = 69.4}");

  tot_filename = (TString) "plotsPU/FinalErrorsPU_" + type + (TString) "_" + method + (TString) ".root"; 
  TFile *f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(finalErrorsPU,"graph");
  f->Close();
  delete f;

  return 0;
}






