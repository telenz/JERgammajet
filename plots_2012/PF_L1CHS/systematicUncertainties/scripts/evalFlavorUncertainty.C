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


int evalFlavorUncertainty(TString definition = "algo"){

  cout<<endl<<endl<<endl<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for Flavor uncertainty is executed! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl;
  gErrorIgnoreLevel = 1001;

  TeresaPlottingStyle::init();
  gStyle -> SetTitleOffset(1.4,"Y");
  gROOT->ForceStyle(); 
  
  const TString method = "RMS99";
  const TString type   = "PFCHS";

  TString name[5];
  name[0]="together";
  name[1]="gluons";
  name[2]="quarks";
  name[3]="uds";
  name[4]="cb";

  
  TString pathName[5];
  for(int i=0; i<5; i++){
    pathName[i]       = (TString) "root_files_FlavorUncertainty/"+name[i]+"_" + definition+ "/";
  }


  TString pathNameTogether = (TString) "root_files_FlavorUncertainty/"+name[0]+"_" + definition+ "/";
  TString pathNameGluons   = (TString) "root_files_FlavorUncertainty/"+name[1]+"_" + definition+ "/";
  TString pathNameQuarks   = (TString) "root_files_FlavorUncertainty/"+name[2]+"_" + definition+ "/";
  TString pathNameUDS      = (TString) "root_files_FlavorUncertainty/"+name[3]+"_" + definition+ "/";
  TString pathNameCB       = (TString) "root_files_FlavorUncertainty/"+name[4]+"_" + definition+ "/";

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Evaluation of the flavor composition
  
  TH2D *charm2D, *bottom2D, *gluon2D, *lightQuark2D;
  TH1D *charm1D[nEtaBins], *bottom1D[nEtaBins], *gluon1D[nEtaBins], *lightQuark1D[nEtaBins];
  double allQuarks[nEtaBins][nPtBins] = {{0.}};
  double togetherWoUndefined[nEtaBins][nPtBins] = {{0.}};
  double allQuarksComp[nEtaBins][nPtBins] = {{0.}};
  double charm[nEtaBins][nPtBins] = {{0.}};
  double bottom[nEtaBins][nPtBins] = {{0.}};
  double lightQuark[nEtaBins][nPtBins] = {{0.}};
  double gluon[nEtaBins][nPtBins] = {{0.}};
  TString fileName;  
  TFile *fileFlavor;  

  double mixture[nEtaBins][nPtBins]    ={{0.}};
  double mixtureUp[nEtaBins][nPtBins]  ={{0.}};
  double mixtureLow[nEtaBins][nPtBins] ={{0.}};

  double correlationQuarks[nEtaBins][nPtBins] = {{0.}};
  double correlationGluons[nEtaBins][nPtBins] = {{0.}};
  
  for(int i=0; i<nEtaBins; i++){
  
    // Read Files
    fileName =  pathName[0] + "hPhotonPtJetEtaFlavorBottom_" + definition + "_PFCHS_mc.root";
    fileFlavor     =  new TFile(fileName);
    bottom2D =  (TH2D*) gDirectory->Get("histo");
    bottom2D -> SetDirectory(0);
    bottom1D[i]  = (TH1D*) bottom2D->ProjectionX("bottom0",bottom2D->GetYaxis()->FindBin(etaBins[i]),bottom2D->GetYaxis()->FindBin(etaBins[i+1]));
    bottom1D[i] -> SetDirectory(0);
    delete fileFlavor;

    fileName =  pathName[0] + "hPhotonPtJetEtaFlavorCharm_" + definition + "_PFCHS_mc.root";
    fileFlavor     =  new TFile(fileName);
    charm2D  =  (TH2D*) gDirectory->Get("histo");
    charm2D  -> SetDirectory(0);
    charm1D[i] = charm2D->ProjectionX("charm0",charm2D->GetYaxis()->FindBin(etaBins[i]),charm2D->GetYaxis()->FindBin(etaBins[i+1]),"");
    charm1D[i]  -> SetDirectory(0);
    delete fileFlavor;

    fileName =  pathName[0] + "hPhotonPtJetEtaFlavorGluon_" + definition + "_PFCHS_mc.root";
    fileFlavor     =  new TFile(fileName);
    gluon2D  =  (TH2D*) gDirectory->Get("histo");
    gluon2D  -> SetDirectory(0);
    gluon1D[i]    = gluon2D->ProjectionX("gluon0",gluon2D->GetYaxis()->FindBin(etaBins[i]),gluon2D->GetYaxis()->FindBin(etaBins[i+1]),"");
    gluon1D[i] -> SetDirectory(0); 
    delete fileFlavor;

    fileName      =  pathName[0] + "hPhotonPtJetEtaFlavorLightQuarks_" + definition + "_PFCHS_mc.root";
    fileFlavor    =  new TFile(fileName);
    lightQuark2D  =  (TH2D*) gDirectory->Get("histo");
    lightQuark2D  -> SetDirectory(0);
    lightQuark1D[i] = lightQuark2D->ProjectionX("lightQuark0",lightQuark2D->GetYaxis()->FindBin(etaBins[i]),lightQuark2D->GetYaxis()->FindBin(etaBins[i+1]),"");
    lightQuark1D[i] -> SetDirectory(0);
    
    
    for(int j=0; j<nPtBins; j++){
      
      for(int bin=lightQuark1D[i]->GetXaxis()->FindBin(ptBins[j]); bin<lightQuark1D[i]->GetXaxis()->FindBin(ptBins[j+1]); bin++){
	lightQuark[i][j] += lightQuark1D[i]->GetBinContent(bin);
      }
      for(int bin=bottom1D[i]->GetXaxis()->FindBin(ptBins[j]); bin<bottom1D[i]->GetXaxis()->FindBin(ptBins[j+1]); bin++){
	bottom[i][j] += bottom1D[i]->GetBinContent(bin);
      }
      for(int bin=charm1D[i]->GetXaxis()->FindBin(ptBins[j]); bin<=charm1D[i]->GetXaxis()->FindBin(ptBins[j+1]); bin++){
	charm[i][j] += charm1D[i]->GetBinContent(bin);
      }
      
      for(int bin=gluon1D[i]->GetXaxis()->FindBin(ptBins[j]); bin<gluon1D[i]->GetXaxis()->FindBin(ptBins[j+1]); bin++){
	gluon[i][j] += gluon1D[i]->GetBinContent(bin);
      }

      allQuarks[i][j] = bottom[i][j];
      allQuarks[i][j] += charm[i][j];
      allQuarks[i][j] += lightQuark[i][j];

      togetherWoUndefined[i][j]  = allQuarks[i][j];
      togetherWoUndefined[i][j] += gluon[i][j];

      correlationQuarks[i][j] = sqrt(allQuarks[i][j]/togetherWoUndefined[i][j]);
      correlationGluons[i][j] = sqrt(gluon[i][j]/togetherWoUndefined[i][j]);
    
      allQuarksComp[i][j]     = allQuarks[i][j]/togetherWoUndefined[i][j];
   
      mixture[i][j]    = allQuarksComp[i][j];
      mixtureUp[i][j]  = allQuarksComp[i][j] + 0.10;
      mixtureLow[i][j] = allQuarksComp[i][j] - 0.10;

    }
    //cout<<"--------------------"<<endl;
  }

  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TString etaString = "Uncertainty of Flavor Composition";
  
  TString tot_filename, AuxString, fitName;;
  TGraphErrors* graph[5]; 
  double finalErrorsUpX[nEtaBins][nPtBins]   = {{0}};
  double finalErrorsUpEX[nEtaBins][nPtBins]  = {{0}};
  double finalErrorsLowX[nEtaBins][nPtBins]  = {{0}};
  double finalErrorsLowEX[nEtaBins][nPtBins] = {{0}};
  double finalErrorsUpY[nEtaBins][nPtBins]   = {{0}};
  double finalErrorsUpEY[nEtaBins][nPtBins]  = {{0}};
  double finalErrorsLowY[nEtaBins][nPtBins]  = {{0}};
  double finalErrorsLowEY[nEtaBins][nPtBins] = {{0}};
  int nCount[nEtaBins] = {0};
  TString rootFile[5];
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  for(int eta = 0; eta<nEtaBins; eta++){

    for(int i=0; i<5; i++){
    rootFile[i]  = pathName[i] + (TString) "Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_intrinsic_" + type + (TString) "_mc_" + method + (TString) ".root";
    }

    TMultiGraph* mg = new TMultiGraph();
    etaString = "Uncertainty on Flavor Composition";
    //mg -> SetTitle(etaString);
    
    TCanvas *c = new TCanvas("c",etaString,200,10,800,800);
    c -> cd();
    
    double maximum = 0;
     
    TString newName;

    for(int i =0; i<5; i++){
      newName = (TString) "Graph_" +name[i];
      graph[i] = readTGraphErrors(rootFile[i],"Graph",newName);    
      if(graph[i]->GetYaxis()->GetXmax() > maximum) maximum = graph[i]->GetYaxis()->GetXmax();
    }

    int colorTogether = kGray+3;
    int colorGluons   = kBlue;
    int colorQuarks   = kGreen+1;
    int colorUDS      = kOrange;
    int colorCB       = kRed;

    
    graph[0] -> SetMarkerStyle(20);
    graph[1] -> SetMarkerStyle(22);
    graph[2] -> SetMarkerStyle(30);
    graph[3] -> SetMarkerStyle(33);
    graph[4] -> SetMarkerStyle(21);

    graph[0] -> SetMarkerSize(1.1);
    graph[1] -> SetMarkerSize(1.3);
    graph[2] -> SetMarkerSize(1.1);
    graph[3] -> SetMarkerSize(1.5);
    graph[4] -> SetMarkerSize(1.0);

    graph[0] -> SetMarkerColor(colorTogether);
    graph[1] -> SetMarkerColor(colorGluons);
    graph[2] -> SetMarkerColor(colorQuarks);
    graph[3] -> SetMarkerColor(colorUDS);
    graph[4] -> SetMarkerColor(colorCB);

    graph[0] -> SetLineColor(colorTogether);
    graph[1] -> SetLineColor(colorGluons);
    graph[2] -> SetLineColor(colorQuarks);
    graph[3] -> SetLineColor(colorUDS);
    graph[4] -> SetLineColor(colorCB);


    TF1* f[5];
   
    for(int i=0; i<5; i++){
      f[i]=graph[i] -> GetFunction("fResolution");
      newName = (TString) "f_" +name[i];
      f[i]->SetName(newName);
    }

    f[0] -> SetLineColor(colorTogether);
    f[1] -> SetLineColor(colorGluons);
    f[2] -> SetLineColor(colorQuarks);
    f[3] -> SetLineColor(colorUDS);
    f[4] -> SetLineColor(colorCB);

    f[0] -> SetFillColor(colorTogether);
    f[1] -> SetFillColor(colorGluons);
    f[2] -> SetFillColor(colorQuarks);
    f[3] -> SetFillColor(colorUDS);
    f[4] -> SetFillColor(colorCB);
    
    graph[0] -> GetFunction("f_together") -> SetLineColor(colorTogether);
    graph[1] -> GetFunction("f_gluons")   -> SetLineColor(colorGluons);
    graph[2] -> GetFunction("f_quarks")   -> SetLineColor(colorQuarks);
    graph[3] -> GetFunction("f_uds")      -> SetLineColor(colorUDS);
    graph[4] -> GetFunction("f_cb")       -> SetLineColor(colorCB);
    


    TLegend *legend  = new TLegend(0.55,0.65,0.9,0.9);
    legend -> SetTextSize(0.045);
    legend -> AddEntry(graph[0],"full sample","pl");
    legend -> AddEntry(graph[1],"gluons","pl");
    legend -> AddEntry(graph[2],"all quarks","pl");
    legend -> AddEntry(graph[3],"uds quarks","pl");
    legend -> AddEntry(graph[4],"cb quarks","pl");
  
    graph[0] -> Draw("AP");
    graph[0] -> SetTitle(etaString);
    graph[0] -> GetYaxis() -> SetTitle("JER");
    graph[0] -> SetMinimum(0.04);
    graph[0] -> SetMaximum(0.14);   
    graph[0] -> GetXaxis() -> SetRangeUser(0.,600.);
    graph[0] -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]"); 

    for(int i=1; i<5; i++) graph[i]->Draw("same P");
  
    
    legend -> Draw("same");
    
    TLatex*  info   = new TLatex();
    info-> SetNDC();
    AuxString = Form("%4.1f < |#eta| < %4.1f",etaBins[eta],etaBins[eta+1]);
    info->DrawLatex(0.65,0.55,AuxString);
    info->DrawLatex(0.2,0.22,"Anti-k_{T} R=0.5");
    info->DrawLatex(0.2,0.16,"PF+CHS");
    
    tot_filename = (TString) "plotsFlavor/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_FlavorUncertainty_" + method ;
    c->Update();
    c -> SaveAs(tot_filename + (TString) ".pdf");
    c -> SaveAs(tot_filename + (TString) ".C");
    
    delete c;
 
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // 2.) Calculate Relative uncertainty for a 10% different mixture with measured resolution

    rootFile[0]  = pathName[0] + (TString) "Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_intrinsic_" + type + (TString) "_mc_" + method + (TString) ".root"; // Together
    rootFile[1]  = pathName[1] + (TString) "Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_intrinsic_" + type + (TString) "_mc_" + method + (TString) ".root"; // Gluons
    rootFile[2]  = pathName[2] + (TString) "Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_intrinsic_" + type + (TString) "_mc_" + method + (TString) ".root"; // Quarks
    
    for(int i =0; i<3; i++){
      graph[i] = readTGraphErrors(rootFile[i],"Graph","Graph"+name[i]);    
      graph[i] -> SetMarkerStyle(20);
      if(graph[i]->GetYaxis()->GetXmax() > maximum) maximum = graph[i]->GetYaxis()->GetXmax();
    }
    
    int nDataGluon = graph[1]->GetN();
    int nDataQuark = graph[2]->GetN();
    int nData = 0;

    if(nDataGluon < nDataQuark) nData = nDataGluon;
    else nData = nDataQuark;



    double *dataY, *dataQuarksY, *dataGluonsY, *dataX, *dataQuarksX, *dataGluonsX, *dataEY, *dataQuarksEY, *dataGluonsEY, *dataEX ;

    dataY        = graph[0] -> GetY();
    dataQuarksY  = graph[2] -> GetY();
    dataGluonsY  = graph[1] -> GetY();
    dataX        = graph[0] -> GetX();
    dataQuarksX  = graph[2] -> GetX();
    dataGluonsX  = graph[1] -> GetX();
    dataEY       = graph[0] -> GetEY();
    dataQuarksEY = graph[2] -> GetEY();
    dataGluonsEY = graph[1] -> GetEY();
    dataEX       = graph[0] -> GetEX();

    double *errorUpY  = new double[nData];
    double *errorLowY = new double[nData];
    double *errorUpEY = new double[nData];
    double *errorLowEY= new double[nData];
    

    int idxGluon   = 0;
    int idxQuark   = 0;
    int idx        = 0;
    int idxMix     = 0;
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

      

      for(int i=0; i<nPtBins; i++){

	if(dataX[idx]>ptBins[i] && dataX[idx]<ptBins[i+1]) idxMix = i;
      }
      
      dataY[idx]            = (1.-mixture[eta][idxMix])*(dataGluonsY[idxGluon])+mixture[eta][idxMix]*(dataQuarksY[idxQuark]);
      errorUpY[countNData]  = ((1.-mixtureUp[eta][idxMix])*(dataGluonsY[idxGluon])+mixtureUp[eta][idxMix]*(dataQuarksY[idxQuark]))/(dataY[idx]) -1.;
      errorLowY[countNData] = ((1.-mixtureLow[eta][idxMix])*dataGluonsY[idxGluon]+mixtureLow[eta][idxMix]*dataQuarksY[idxQuark])/(dataY[idx]) -1.;

      errorUpEY[countNData] = sqrt(pow(1./(dataY[idx]),2)*(pow((1.-mixtureUp[eta][idxMix])*dataGluonsEY[idxGluon],2) + pow(mixtureUp[eta][idxMix]*dataQuarksEY[idxQuark],2)) + pow(abs((1.-mixtureUp[eta][idxMix])*dataGluonsY[idxGluon]+mixtureUp[eta][idxMix]*dataQuarksY[idxQuark])/(pow(dataY[idx],2)),2)*pow(dataEY[idx],2) - ((1.-mixtureUp[eta][idxMix])*correlationGluons[eta][idxMix]*dataGluonsEY[idxGluon] + mixtureUp[eta][idxMix]*correlationQuarks[eta][idxMix]*dataQuarksEY[idxQuark])*2.0 * abs((1.-mixtureUp[eta][idxMix])*dataGluonsY[idxGluon]+mixtureUp[eta][idxMix]*dataQuarksY[idxQuark])/(pow(dataY[idx],3))*dataEY[idx]);
            
      errorLowEY[countNData] = sqrt(pow(1./(dataY[idx]),2)*(pow((1.-mixtureLow[eta][idxMix])*dataGluonsEY[idxGluon],2) + pow(mixtureLow[eta][idxMix]*dataQuarksEY[idxQuark],2)) + pow(abs((1.-mixtureLow[eta][idxMix])*dataGluonsY[idxGluon]+mixtureLow[eta][idxMix]*dataQuarksY[idxQuark])/(pow(dataY[idx],2)),2)*pow(dataEY[idx],2) - ((1.-mixtureLow[eta][idxMix])*correlationGluons[eta][idxMix]*dataGluonsEY[idxGluon] + mixtureLow[eta][idxMix]*correlationQuarks[eta][idxMix]*dataQuarksEY[idxQuark])*2.0 * abs((1.-mixtureLow[eta][idxMix])*dataGluonsY[idxGluon]+mixtureLow[eta][idxMix]*dataQuarksY[idxQuark])/(pow(dataY[idx],3))*dataEY[idx]);

      dataX[countNData]  = dataX[idx];
      dataEX[countNData] = dataEX[idx];
      
      countNData += 1;
      idx        += 1;
      idxQuark   += 1;
      idxGluon   += 1;
    }

    TGraphErrors *plotUp  = new TGraphErrors(nData,dataX,errorUpY,dataEX,errorUpEY);
    TGraphErrors *plotLow = new TGraphErrors(nData,dataX,errorLowY,dataEX,errorLowEY);

    plotUp  -> GetXaxis() -> SetRangeUser(0.,600.);  
    plotLow -> GetXaxis() -> SetRangeUser(0.,600.);  

    TF1 *Line   =  new TF1("Line","pol0",0.,600.);

    Line->SetParameter(0,0);

    TCanvas* c1 = new TCanvas("c1","c1",200,10,800,800);
    c1 -> cd();

    plotUp -> SetMarkerStyle(20);  
    plotUp -> SetMarkerColor(1);  
    plotUp -> SetLineColor(1); 
    
    plotLow -> SetMarkerStyle(24);  
    plotLow -> SetMarkerColor(9);  
    plotLow -> SetLineColor(9);
 
    
    TLegend* legend1  = new TLegend(0.30,0.75,0.9,0.9);
    legend1 -> SetTextSize(0.045);
    
    legend1 -> AddEntry(plotLow,Form("- 10%% quarks/ + 10%% gluons"),"p");
    legend1 -> AddEntry(plotUp ,Form("+ 10%% quarks/ - 10%% gluons"),"p");
 
    mg = new TMultiGraph();
    mg->Add(plotLow);  
    mg->Add(plotUp);

    etaString = "Relative Flavor Uncertainty";
    mg -> SetTitle(etaString);  
  
    mg->Draw("AP");

    mg -> GetYaxis() -> SetTitle("JER_{different mixture} /JER_{full sample} - 1");
    mg -> SetMinimum(-0.1);
    mg -> SetMaximum(0.1); 
    mg -> GetXaxis() -> SetLimits(0.,600.);  
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  
    legend1->Draw("same");

    Line->Draw("same");

    
    TLatex* info1   = new TLatex();
    info1-> SetNDC();
    info1->DrawLatex(0.6,0.7,AuxString);
  
    double *xUp   = plotUp  -> GetX();
    double *xEUp  = plotUp  -> GetEX();
    double *xLow  = plotLow -> GetX();
    double *xELow = plotLow -> GetEX();
    double *yUp   = plotUp  -> GetY();
    double *yEUp  = plotUp  -> GetEY();
    double *yLow  = plotLow -> GetY();
    double *yELow = plotLow -> GetEY();


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
      finalErrorsUpEX[eta][pt] = xEUp[pt];
      finalErrorsUpY[eta][pt]  = yUp[pt];
      finalErrorsUpEY[eta][pt] = yEUp[pt];
      finalErrorsLowX[eta][pt] = xLow[pt];
      finalErrorsLowEX[eta][pt]= xELow[pt];
      finalErrorsLowY[eta][pt] = yLow[pt];
      finalErrorsLowEY[eta][pt]= yELow[pt];

      nCount[eta] += 1;


    }

    tot_filename = (TString) "plotsFlavor/Relative_Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_FlavorUncertainty_" + method + (TString) "_mixture.pdf";
  
    c1 -> SaveAs(tot_filename);

    delete mg;  
    delete c1;
  

  }

  // Save relative uncertainties for every eta bin in another root-file

  TGraphErrors* finalErrorsFlavorUp[nEtaBins];
  TGraphErrors* finalErrorsFlavorLow[nEtaBins]; 


  for(int i=0; i<nEtaBins; i++){


    finalErrorsFlavorUp[i] = new TGraphErrors(nCount[i],finalErrorsUpX[i],finalErrorsUpY[i],finalErrorsUpEX[i],finalErrorsUpEY[i]);
    finalErrorsFlavorUp[i] -> SetMarkerStyle(20);
    finalErrorsFlavorUp[i] -> SetTitle("Final relative Erros (Flavor)");
    finalErrorsFlavorUp[i] -> GetXaxis() -> SetTitle("#eta^{Jet}");
    finalErrorsFlavorUp[i] -> GetYaxis() -> SetTitle("JER_{quarks/gluons} /JER_{full sample}");
    
    finalErrorsFlavorLow[i] = new TGraphErrors(nCount[i],finalErrorsLowX[i],finalErrorsLowY[i],finalErrorsLowEX[i],finalErrorsLowEY[i]);
    finalErrorsFlavorLow[i] -> SetMarkerStyle(20);
    finalErrorsFlavorLow[i] -> SetTitle("Final relative Erros (Flavor)");
    finalErrorsFlavorLow[i] -> GetXaxis() -> SetTitle("#eta^{Jet}");
    finalErrorsFlavorLow[i] -> GetYaxis() -> SetTitle("JER_{quarks/gluons} /JER_{full sample}");

    tot_filename = (TString) "plotsFlavor/RelativeErrorsPtBinnedFlavorUp_" + type + (TString) "_" + (long) (i+1) + "_bin_" + method + (TString) ".root"; 
    TFile *f = new TFile(tot_filename,"RECREATE");
    f -> WriteTObject(finalErrorsFlavorUp[i],"graph");
    f->Close();
    delete f;
    tot_filename = (TString) "plotsFlavor/RelativeErrorsPtBinnedFlavorLow_" + type + (TString) "_" + (long) (i+1) + "_bin_" + method + (TString) ".root"; 
    f = new TFile(tot_filename,"RECREATE");
    f -> WriteTObject(finalErrorsFlavorLow[i],"graph");
    f->Close();
    delete f;
  }

  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Move MC points up and down with the errors and fit again the data/MC_up and data/MC_down ratios -> Get errors for every eta Bin

  TString rootFiles;
  TString JetType = "PFCHS";

  double *ratioEtaBinnedX     = new double[nEtaBins];
  double *ratioEtaBinnedY     = new double[nEtaBins];
  double *ratioEtaBinnedUpY   = new double[nEtaBins];
  double *ratioEtaBinnedLowY  = new double[nEtaBins];
  double *ratioEtaBinnedEX    = new double[nEtaBins];
  double *ratioEtaBinnedEY    = new double[nEtaBins];
  double *ratioEtaBinnedUpEY  = new double[nEtaBins];
  double *ratioEtaBinnedLowEY = new double[nEtaBins];

  for(int iEta = 0; iEta < nEtaBins; iEta++){
    
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
    double *mcUpY   = new double[nData];
    double *mcLowY  = new double[nData];
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

      JERMC    -> GetPoint(idxMC,mcX[0],mcY[0]);
      mcUpY[0]  = mcY[0];
      mcLowY[0] = mcY[0];
      mcEX[0]   = JERMC -> GetErrorX(idxMC);
      mcEY[0]   = JERMC -> GetErrorY(idxMC);


      if(idxMC>nPtBins || idxUp>nPtBins || idxLow>nPtBins) break;

      if(abs(mcX[0]/finalErrorsUpX[iEta][idxUp]-1.)>0.05){
	idxMC += 1;
	continue;
      }
      else{
	mcUpY[0] = mcY[0]*(1.+finalErrorsUpY[iEta][idxUp]);
      }

      if(abs(mcX[0]/finalErrorsLowX[iEta][idxLow]-1.)>0.05){
	idxMC += 1;
	continue;
      }
      else{
	mcLowY[0] = mcY[0]*(1.+finalErrorsLowY[iEta][idxLow]);
      }

      if(abs(dataX[idxData]/mcX[0] - 1.) > 0.05){

	if(dataX[idxData]<mcX[0]){
	  idxData += 1;
	  nData -= 1;
	  continue; 
	}
	else{
	  idxMC  += 1;
	  idxUp  += 1;
	  idxLow += 1;
	  continue;
	}
      }
          
      ratioX[idx]  = 1./2.*(dataX[idxData] + mcX[0]);
      ratioY[idx]  = dataY[idxData]/mcY[0];
      ratioUpY[idx]   = dataY[idxData]/mcLowY[0];
      ratioLowY[idx]  = dataY[idxData]/mcUpY[0];
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
    Ratio -> GetXaxis() -> SetTitle("Photon pT");
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


    
    TString filename = (TString) "plotsFlavor/Ratios_upper_lower_Error_for_" + (long) (iEta+1) + (TString) "_eta_bin_" + type + (TString) "_" + method + (TString) ".pdf";
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

  TGraphErrors* ratioEtaBinned    = new TGraphErrors(nEtaBins,ratioEtaBinnedX,ratioEtaBinnedY,ratioEtaBinnedEX,ratioEtaBinnedEY);
  TGraphErrors* ratioEtaBinnedUp  = new TGraphErrors(nEtaBins,ratioEtaBinnedX,ratioEtaBinnedUpY,ratioEtaBinnedEX,ratioEtaBinnedUpEY);
  TGraphErrors* ratioEtaBinnedLow = new TGraphErrors(nEtaBins,ratioEtaBinnedX,ratioEtaBinnedLowY,ratioEtaBinnedEX,ratioEtaBinnedLowEY);

  TString filename = (TString) "plotsFlavor/RatioEtaBinned_" + type + (TString) "_" + method + (TString) ".root";
  TFile *f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinned,"Graph");
  f->Close();
  delete f;

  filename = (TString) "plotsFlavor/FinalEtaBinnedErrorsFlavorUp_" + type + (TString) "_"  + method + (TString) ".root"; 
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinnedUp,"graph");
  f->Close();
  delete f;

  filename = (TString) "plotsFlavor/FinalEtaBinnedErrorsFlavorLow_" + type + (TString) "_"  + method + (TString) ".root"; 
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinnedLow,"graph");
  f->Close();
  delete f;

  return 0;
 
}
