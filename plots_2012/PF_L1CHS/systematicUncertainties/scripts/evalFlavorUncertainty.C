// Script to evaluate Flavor Uncertainty

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


int evalFlavorUncertainty(TString definition){

  TeresaPlottingStyle::init();

  const int nEta = 4;
  
  const TString method = "RMS99";
  const TString type   = "PFCHS";

  //TString definition = "algo";
  //TString definition = "phys";

  TString pathName = "root_files_WithoutTriggerWithPUWeightEq1/";
  TString pathNameGluons = (TString) "root_files_FlavorUncertainty/gluons_" + definition+ "/";
  TString pathNameQuarks = (TString) "root_files_FlavorUncertainty/quarks_" + definition+ "/";


  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Evaluation of the flavor composition
  
  TH2D *charm2D, *bottom2D, *gluon2D, *lightQuark2D;
  TH1D *charm[2], *bottom[2], *gluon[2], *lightQuark[2],  *allQuarks[2], *togetherWoUndefined[2], *allQuarksComp[2], *gluonWoUndefinedComp[2];
  TString fileName;  
  TFile *fileFlavor;  
   
  // Read Files
  fileName =  pathName + "hPhotonPtJetEtaFlavorBottom_" + definition + "_PFCHS_mc.root";
  fileFlavor     =  new TFile(fileName);
  bottom2D =  (TH2D*) gDirectory->Get("histo");
  bottom2D -> SetDirectory(0);
  bottom[0]  = (TH1D*) bottom2D->ProjectionX("bottom0",0,bottom2D->GetYaxis()->FindBin(1.3));
  bottom[1]  = (TH1D*) bottom2D->ProjectionX("bottom1",bottom2D->GetYaxis()->FindBin(1.3),-1);
  for(int i= 0; i<2; i++) bottom[i] -> Rebin(20);
  for(int i= 0; i<2; i++) bottom[i] -> SetDirectory(0);
  delete fileFlavor;

  fileName =  pathName + "hPhotonPtJetEtaFlavorCharm_" + definition + "_PFCHS_mc.root";
  fileFlavor     =  new TFile(fileName);
  charm2D  =  (TH2D*) gDirectory->Get("histo");
  charm2D  -> SetDirectory(0);
  charm[0] = charm2D->ProjectionX("charm0",0,charm2D->GetYaxis()->FindBin(1.3),"");  
  charm[1] = charm2D->ProjectionX("charm1",charm2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++) charm[i]  ->Rebin(20);
  for(int i= 0; i<2; i++) charm[i]  -> SetDirectory(0);
  delete fileFlavor;

  fileName =  pathName + "hPhotonPtJetEtaFlavorGluon_" + definition + "_PFCHS_mc.root";
  fileFlavor     =  new TFile(fileName);
  gluon2D  =  (TH2D*) gDirectory->Get("histo");
  gluon2D  -> SetDirectory(0);
  gluon[0]    = gluon2D->ProjectionX("gluon0",0,gluon2D->GetYaxis()->FindBin(1.3),"");
  gluon[1]    = gluon2D->ProjectionX("gluon1",gluon2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++)  gluon[i]    -> Rebin(20);
  for(int i= 0; i<2; i++)  gluon[i]      -> SetDirectory(0); 
  delete fileFlavor;
  
  fileName     =  pathName + "hPhotonPtJetEtaFlavorLightQuarks_" + definition + "_PFCHS_mc.root";
  fileFlavor         =  new TFile(fileName);
  lightQuark2D =  (TH2D*) gDirectory->Get("histo");
  lightQuark2D -> SetDirectory(0);
  lightQuark[0] = lightQuark2D->ProjectionX("lightQuark0",0,lightQuark2D->GetYaxis()->FindBin(1.3),"");
  lightQuark[1] = lightQuark2D->ProjectionX("lightQuark1",lightQuark2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++)   lightQuark[i] -> Rebin(20);
  for(int i= 0; i<2; i++)   lightQuark[i] -> SetDirectory(0);

     
  for(int i= 0; i<1; i++){

    allQuarks[i] = (TH1D*) bottom[i]->Clone(Form("allQuarks%i",i));
    allQuarks[i] -> Add(charm[i]);
    allQuarks[i] -> Add(lightQuark[i]);

    togetherWoUndefined[i] = (TH1D*) bottom[i]->Clone(Form("together%i",i));
    togetherWoUndefined[i] -> Add(charm[i]);
    togetherWoUndefined[i] -> Add(gluon[i]);
    togetherWoUndefined[i] -> Add(lightQuark[i]);

    togetherWoUndefined[i]-> Rebin(togetherWoUndefined[i]->GetNbinsX());
    allQuarks[i]          -> Rebin(allQuarks[i]->GetNbinsX());
    gluon[i]              -> Rebin(gluon[i]->GetNbinsX());
    
    allQuarksComp[i]        = (TH1D*) allQuarks[i]->Clone(Form("allQuarks%i",i+3));
    allQuarksComp[i]        -> Divide(togetherWoUndefined[i]);    
    gluonWoUndefinedComp[i] = (TH1D*) gluon[i]->Clone(Form("gluon%i",i+3));
    gluonWoUndefinedComp[i] -> Divide(togetherWoUndefined[i]);

    cout<<endl<<endl<<"Gluon Flavor Fraction for systematics  = "<<gluonWoUndefinedComp[i]->GetBinContent(1)<<endl;
    cout<<"Quark Flavor Fraction for systematics  = "<<allQuarksComp[i]->GetBinContent(1)<<endl;

  }
  double mixture    = allQuarksComp[0]->GetBinContent(1);
  double mixtureUp  = allQuarksComp[0]->GetBinContent(1) + 0.10;
  double mixtureLow = allQuarksComp[0]->GetBinContent(1) - 0.10;   

  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  TString etaString = "Uncertainty on Flavor Composition";
  
  TString tot_filename, AuxString, fitName;;
  TGraphErrors* graph[3]; 
  double finalErrorsUp[nEta]   = {0};
  double finalErrorsUpE[nEta]  = {0};
  double finalErrorsLow[nEta]  = {0};
  double finalErrorsLowE[nEta] = {0};


  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Calculate Correlation with the help of the control plot hPhoton1Pt_PFCHS_mc.root
  TString rootFile[3]; 
  rootFile[0]  =  pathName + (TString) "hPhoton1Pt_" + type + (TString) "_mc.root";
  rootFile[1]  =  pathNameQuarks + "hPhoton1Pt_" + type + (TString) "_mc.root";
  rootFile[2]  =  pathNameGluons + "hPhoton1Pt_" + type + (TString) "_mc.root"; 

  TFile* file[3];
  TH1D* histo[3];
  
  for(int i =0; i<3; i++){
    file[i] = TFile::Open(rootFile[i]);
    file[i]->GetObject("histo",histo[i]);    
  }
   
  const double correlationQuarks = sqrt(histo[1]->Integral()/histo[0]->Integral());
  const double correlationGluons = sqrt(histo[2]->Integral()/histo[0]->Integral());

  cout<<endl<<"Correlation between Quark and full Sample = "<<correlationQuarks<<endl;
  cout<<      "Correlation between Gluon and full Sample = "<<correlationGluons<<endl<<endl;
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  for(int eta = 1; eta<nEta+1; eta++){
    
    rootFile[0]  = pathName + (TString) "Resolution_for_" + (long) eta + (TString) "_eta_bin_" + type + (TString) "_mc_" + method + (TString) ".root";
    rootFile[1]  = pathNameQuarks + (TString) "Resolution_for_" + (long) eta + (TString) "_eta_bin_" + type + (TString) "_mc_" + method + (TString) ".root";
    rootFile[2]  = pathNameGluons + (TString) "Resolution_for_" + (long) eta + (TString) "_eta_bin_" + type + (TString) "_mc_" + method + (TString) ".root";    
  
    TMultiGraph* mg = new TMultiGraph();
    etaString = "Uncertainty on Flavor Composition";
    mg -> SetTitle(etaString);
    
    TCanvas *c = new TCanvas("c",etaString,200,10,800,800);
    c -> cd();
    
    double maximum = 0;
     
    for(int i =0; i<3; i++){
      graph[i] = readTGraphErrors(rootFile[i],"Graph","Graph");    
      graph[i] -> SetMarkerStyle(20);
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
  
 
    TLegend *legend  = new TLegend(0.35,0.75,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.042);
    legend -> AddEntry(graph[0],"Full sample","l");
    legend -> AddEntry(graph[1],"Quarks","l");
    legend -> AddEntry(graph[2],"Gluons","l");
  
    for(int i=0; i<3; i++) mg->Add(graph[i]);  
  
    mg->Draw("AP");
    
    mg -> GetYaxis() -> SetTitle("JER");
    mg -> SetMinimum(0.004);
    mg -> SetMaximum(0.12);   
    mg -> GetXaxis() -> SetLimits(0,600);
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]"); 
    legend -> Draw("same");

    TLatex*  info   = new TLatex();
    info-> SetNDC();
    AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta-1],etaBins[eta]);
    info->DrawLatex(0.6,0.7,AuxString);

    tot_filename = (TString) "plotsFlavor/Resolution_for_" + (long) eta + (TString) "_eta_bin_FlavorUncertainty_" + method + (TString) ".pdf";
    c -> SaveAs(tot_filename);
 
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // 2.) Calculate Relative uncertainty for a 10% different mixture
    
    int nDataQuark = graph[1]->GetN();
    int nDataGluon = graph[2]->GetN();
    int nData = 0;
    if(nDataGluon < nDataQuark) nData = nDataGluon;
    else nData = nDataQuark;


    double *dataY, *dataQuarksY, *dataGluonsY, *dataX, *dataQuarksX, *dataGluonsX, *dataEY, *dataQuarksEY, *dataGluonsEY, *dataEX ;

    dataY        = graph[0] -> GetY();
    dataQuarksY  = graph[1] -> GetY();
    dataGluonsY  = graph[2] -> GetY();
    dataX        = graph[0] -> GetX();
    dataQuarksX  = graph[1] -> GetX();
    dataGluonsX  = graph[2] -> GetX();
    dataEY       = graph[0] -> GetEY();
    dataQuarksEY = graph[1] -> GetEY();
    dataGluonsEY = graph[2] -> GetEY();
    dataEX       = graph[0] -> GetEX();

    double *errorUpY  = new double[nData];
    double *errorLowY = new double[nData];
    double *errorUpEY = new double[nData];
    double *errorLowEY= new double[nData];
    

    int idxGluon   = 0;
    int idxQuark   = 0;
    int idx        = 0;
    int countNData = 0;

    while(countNData < nData){


      if(TMath::Abs(dataX[idx]/dataGluonsX[idxGluon] - 1.) > 0.05){
	if(dataX[idx] > dataGluonsX[idxGluon]) idxGluon += 1;
	else idx += 1;
	continue;
      }
      if(TMath::Abs(dataX[idx]/dataQuarksX[idxQuark] - 1.) > 0.05){
	if(dataX[idx] > dataQuarksX[idxQuark]) idxQuark += 1;
	else idx += 1;
	continue;
      }

      dataY[idx] = (1.-mixture)*(dataGluonsY[idxGluon])+mixture*(dataQuarksY[idxQuark]);
      //errorUpY[countNData]  = sqrt((1.-mixtureUp)*pow(dataGluonsY[idxGluon],2)+mixtureUp*pow(dataQuarksY[idxQuark],2))/(dataY[idx]) -1.;
      errorUpY[countNData]  = ((1.-mixtureUp)*(dataGluonsY[idxGluon])+mixtureUp*(dataQuarksY[idxQuark]))/(dataY[idx]) -1.;
      errorUpY[countNData]  = ((1.-mixtureUp)*(dataGluonsY[idxGluon])+mixtureUp*(dataQuarksY[idxQuark]))/dataY[idx] -1.;
      errorLowY[countNData] = sqrt((1.-mixtureLow)*pow(dataGluonsY[idxGluon],2)+mixtureLow*pow(dataQuarksY[idxQuark],2))/(dataY[idx]) -1.;

      errorUpEY[countNData] = sqrt(pow(1./(dataY[idx]),2)*(pow((1.-mixtureUp)*dataGluonsEY[idxGluon],2) + pow(mixtureUp*dataQuarksEY[idxQuark],2)) + pow(abs((1.-mixtureUp)*dataGluonsY[idxGluon]-mixtureUp*dataQuarksY[idxQuark])/(pow(dataY[idx],2)),2)*pow(dataEY[idx],2) - ((1.-mixtureUp)*correlationGluons*dataGluonsEY[idxGluon] + mixtureUp*correlationQuarks*dataQuarksEY[idxQuark])*2.0 * abs(dataGluonsY[idxGluon]-dataQuarksY[idxQuark])/(2.*pow(dataY[idx],3))*dataEY[idx]);
      
      errorLowEY[countNData] = sqrt(pow(1./(dataY[idx]),2)*(pow((1.-mixtureLow)*dataGluonsEY[idxGluon],2) + pow(mixtureLow*dataQuarksEY[idxQuark],2)) + pow(abs((1.-mixtureLow)*dataGluonsY[idxGluon]-mixtureLow*dataQuarksY[idxQuark])/(pow(dataY[idx],2)),2)*pow(dataEY[idx],2) - ((1.-mixtureLow)*correlationGluons*dataGluonsEY[idxGluon] + mixtureLow*correlationQuarks*dataQuarksEY[idxQuark])*2.0 * abs(dataGluonsY[idxGluon]-dataQuarksY[idxQuark])/(2.*pow(dataY[idx],3))*dataEY[idx]);

      dataX[countNData]  = dataX[idx];
      dataEX[countNData] = dataEX[idx];

      
      countNData += 1;
      idx        += 1;
      idxQuark   += 1;
      idxGluon   += 1;
    }

    TGraphErrors *plotUp  = new TGraphErrors(nData,dataX,errorUpY,dataEX,errorUpEY);
    TGraphErrors *plotLow = new TGraphErrors(nData,dataX,errorLowY,dataEX,errorLowEY);


    TF1 *fitLineUp   =  new TF1("fitLineUp","pol0",0,600);
    fitLineUp        -> SetLineColor(1);
    plotUp -> Fit("fitLineUp","QR");

    TF1 *fitLineLow   =  new TF1("fitLineLow","pol0",0,600);
    fitLineLow        -> SetLineColor(9);
    plotLow -> Fit("fitLineLow","QR");

    TCanvas* c1 = new TCanvas("c1","c1",200,10,800,800);
    c1 -> cd();

    plotUp -> SetMarkerStyle(20);  
    plotUp -> SetMarkerColor(1);  
    plotUp -> SetLineColor(1); 
    
    plotLow -> SetMarkerStyle(24);  
    plotLow -> SetMarkerColor(9);  
    plotLow -> SetLineColor(9);
 
    
    TLegend* legend1  = new TLegend(0.35,0.75,0.9,0.9);
    legend1 -> SetFillColor(0);
    legend1 -> SetTextSize(0.033);
   
    legend1 -> AddEntry(plotLow,Form("%4.2f #upoint Quarks, %4.2f #upoint Gluons", mixtureLow, 1.-mixtureLow),"p");
    legend1 -> AddEntry(plotUp,Form("%4.2f #upoint Quarks, %4.2f #upoint Gluons", mixtureUp, 1.-mixtureUp),"p");
 
    mg = new TMultiGraph();
    mg->Add(plotLow);  
    mg->Add(plotUp);

    etaString = "Relative Flavor Uncertainty";
    mg -> SetTitle(etaString);  
  
    mg->Draw("AP");


    mg -> GetYaxis() -> SetTitle("JER_{different mixture} /JER_{full sample} ");
    mg -> SetMinimum(-0.2);
    mg -> SetMaximum(0.2); 
    mg -> GetXaxis() -> SetLimits(0,600);  
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  
    legend1->Draw("same");

    
    TLatex* info1   = new TLatex();
    info1->SetTextFont(132);
    info1-> SetNDC();
    info1->SetTextSize(0.040);
    info1->DrawLatex(0.6,0.7,AuxString);
  
    AuxString = Form("Chi2/ndof = %4.1f/%i",fitLineUp->GetChisquare(),fitLineUp->GetNDF());
    info1->DrawLatex(0.6,0.2,AuxString);

    AuxString = Form("f = %4.3f #pm %4.3f",fitLineUp->GetParameter(0),fitLineUp->GetParError(0));
    info1->DrawLatex(0.6,0.25,AuxString);

    AuxString = Form("Chi2/ndof = %4.1f/%i",fitLineLow->GetChisquare(),fitLineLow->GetNDF());
    info1->DrawLatex(0.2,0.2,AuxString);

    AuxString = Form("f = %4.3f #pm %4.3f",fitLineLow->GetParameter(0),fitLineLow->GetParError(0));
    info1->DrawLatex(0.2,0.25,AuxString);

    finalErrorsUp[eta-1]  = fitLineUp -> GetParameter(0);
    finalErrorsUpE[eta-1] = fitLineUp -> GetParError(0);
    finalErrorsLow[eta-1] = fitLineLow -> GetParameter(0);
    finalErrorsLowE[eta-1]= fitLineLow -> GetParError(0);

    tot_filename = (TString) "plotsFlavor/Relative_Resolution_for_" + (long) eta + (TString) "_eta_bin_FlavorUncertainty_" + method + (TString) "_mixture.pdf";
  
    c1 -> SaveAs(tot_filename);

    delete mg;  
    delete c1;
  

  }

  // Save relative uncertainties for every eta bin in another root-file
  double eta[nEta] = {1.};
  double etaError[nEta] = {0.};
  for(int i =0; i<nEta-1; i++) eta[i+1] = eta[i]+1.;

  
  TGraphErrors* finalErrorsFlavorUp = new TGraphErrors(nEta,eta,finalErrorsUp,etaError,finalErrorsUpE);
  finalErrorsFlavorUp -> SetMarkerStyle(20);
  finalErrorsFlavorUp -> SetTitle("Final relative Erros (Flavor)");
  finalErrorsFlavorUp -> GetXaxis() -> SetTitle("#eta^{Jet}");
  finalErrorsFlavorUp -> GetYaxis() -> SetTitle("JER_{quarks/gluons} /JER_{full sample}");

  TGraphErrors* finalErrorsFlavorLow = new TGraphErrors(nEta,eta,finalErrorsLow,etaError,finalErrorsLowE);
  finalErrorsFlavorLow -> SetMarkerStyle(20);
  finalErrorsFlavorLow -> SetTitle("Final relative Erros (Flavor)");
  finalErrorsFlavorLow -> GetXaxis() -> SetTitle("#eta^{Jet}");
  finalErrorsFlavorLow -> GetYaxis() -> SetTitle("JER_{quarks/gluons} /JER_{full sample}");


  tot_filename = (TString) "plotsFlavor/FinalErrorsFlavorUp_" + type + (TString) "_" + method + (TString) ".root"; 
  TFile *f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(finalErrorsFlavorUp,"graph");
  f->Close();
  delete f;
  tot_filename = (TString) "plotsFlavor/FinalErrorsFlavorLow_" + type + (TString) "_" + method + (TString) ".root"; 
  f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(finalErrorsFlavorLow,"graph");
  f->Close();
  delete f;

  return 0;
 
}
