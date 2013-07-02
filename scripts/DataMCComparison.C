// ------------------------------------------------------------------------------------------------
// -------  Script to compare Jet Energy Resolutions and Scales between Data and MC (T.L.)  ------- 
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <iomanip>
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
#include "TVirtualFitter.h"
#include "../CODE/myDeclarations.h"
#include "/afs/naf.desy.de/user/t/telenz/comparison/tdrstyle_mod.C"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);



int plotDataMCComparison(int eta){


  const TString method = "gaus";


  cout<<endl<<"Method is: "<<method<<endl<<endl<<endl;
 
  TString rootFileMC   = (TString) "../plots_2012/PF_L1CHS/mc/root_files/Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_mc_" + method + (TString) ".root";
  TString rootFileData = (TString) "../plots_2012/PF_L1CHS/data/root_files/Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_data_" + method + (TString) ".root";

  cout<<endl<<"Input files are:"<<endl;
  cout<<"Data :"<<rootFileData<<endl;
  cout<<"MC   :"<<rootFileMC<<endl<<endl;

  TString etaString = Form("Jet Energy Resolution");
  
    
  TGraphErrors* JERMC = readTGraphErrors(rootFileMC,"Graph;1","Graph;1");  
  JERMC -> SetMarkerColor(3);
  JERMC -> SetMarkerStyle(20);
  JERMC -> SetLineColor(3);    
  TGraphErrors* JERData = readTGraphErrors(rootFileData,"Graph;1","Graph;1");
  JERData -> SetMarkerColor(4);
  JERData -> SetMarkerStyle(20);
  JERData -> SetLineColor(4);
  
  TCanvas *c = new TCanvas("c",etaString,800,10,700,700);
  c -> cd();
  c -> SetBottomMargin(0.14);
  c -> SetLeftMargin(0.14);
  
  TMultiGraph* mg = new TMultiGraph();
  mg -> SetTitle(etaString);

  TString fitName = "fResolution";
  JERMC   -> GetFunction(fitName) -> SetLineColor(3);
  JERData -> GetFunction(fitName) -> SetLineColor(4);
  

  TLegend *legend  = new TLegend(0.6,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.045);
  legend -> AddEntry(JERData,"Data","l");
  legend -> AddEntry(JERMC,"MC","l");
  
  mg -> Add(JERMC);
  mg -> Add(JERData);
  mg -> Draw("AP");
  
  mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma}");
  mg -> GetYaxis() -> SetTitleOffset(1.3); 
  mg -> GetXaxis() -> SetTitleSize(0.05);
  mg -> GetYaxis() -> SetTitleSize(0.05);
  
  
  mg -> GetYaxis() -> SetTitle("JER");
  mg -> SetMinimum(0.0);
  mg -> SetMaximum(0.2);   

  TLatex*  info   = new TLatex();
  info->SetTextFont(132);
  info-> SetNDC();
  info->SetTextSize(0.048);
  TString AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta-1],etaBins[eta]);
  info->DrawLatex(0.6,0.7,AuxString);
  
  legend -> Draw("same");
 
  TString pdfFile = (TString) "plots/Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_data_mc_comparison_" + method + (TString) ".pdf";
  c -> SaveAs(pdfFile);

  // -------------------  Ratio Plot for Resolution  -------------------
  int nData    = JERData->GetN();
  cout<<"nData = "<<nData<<endl;

  double *dataX  = JERData->GetX();
  double *dataY  = JERData->GetY();
  double *dataEX = JERData->GetEX();
  double *dataEY = JERData->GetEY();
    
  double *mcX  = JERMC->GetX();
  double *mcY  = JERMC->GetY();
  double *mcEX = JERMC->GetEX();
  double *mcEY = JERMC->GetEY();
    
  double *ratioX  = new double[nData];
  double *ratioY  = new double[nData];
  double *ratioEX = new double[nData];
  double *ratioEY = new double[nData];

 
  int idxMC = 0;
  for(int i=0; i<nData; i++){


    if(TMath::Abs(dataX[i]/mcX[idxMC] - 1.) > 0.05){
      i -= 1;
      idxMC += 1;
      continue;
    }
     
    ratioX[i]  = 1./2.*(dataX[i] + mcX[idxMC]);
    ratioY[i]  = dataY[i]/mcY[idxMC];
    ratioEX[i] = 1./2.*TMath::Sqrt(TMath::Power(dataEX[i],2)+TMath::Power(mcEX[idxMC],2));
    ratioEY[i] = TMath::Sqrt(TMath::Power((1./mcY[idxMC]),2)*TMath::Power(dataEY[i],2)+TMath::Power((dataY[i]/(TMath::Power(mcY[idxMC],2))),2)*TMath::Power(mcEY[idxMC],2));
  }
 
  //TGraphErrors *Ratio = new TGraphErrors(10,axcomp,aycomp,axcompError,aycompError);
  TGraphErrors *Ratio = new TGraphErrors(JERData->GetN(),ratioX,ratioY,ratioEX,ratioEY);
  Ratio -> GetXaxis()->SetLimits(0,600);

  TString ratioName = Form("Data/MC Ratio");

  Ratio -> SetTitle(ratioName); 
  Ratio -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  Ratio -> GetXaxis() -> SetTitleOffset(1.1); 
  //Ratio -> GetYaxis() -> SetTitle("Ratio of JER (MC_{smeared}/MC)");
  Ratio -> GetYaxis() -> SetTitle("Ratio of JER (Data/MC)");
  Ratio -> GetYaxis() -> SetTitleOffset(1.2);   

  TF1* f1 = new TF1("name","pol0",0,600);   
  Ratio -> Fit("name","Q");
  cout<<"ChiSquare = "<<f1 -> GetChisquare()<<endl;
  legend  = 0;
  legend = new TLegend(0.55,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.05);
  char legname[100];
  sprintf(legname," %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0));
  legend -> SetHeader(legname);
  TCanvas *c11 = new TCanvas("c11",ratioName,800,10,700,700);
  c11 -> cd();
  c11 -> SetBottomMargin(0.14);
  c11 -> SetLeftMargin(0.14);
  Ratio -> SetMinimum(0.5);
  Ratio -> SetMaximum(2.0);
  

  // Draw info boxes
  Ratio  -> Draw("AP"); 
  Ratio  -> GetXaxis() -> SetTitleSize(0.05);
  Ratio  -> GetYaxis() -> SetTitleSize(0.05);
  legend -> Draw("same");

  TLatex*  info1  = new TLatex();  
  info1 ->SetTextFont(132);
  info1 -> SetNDC();
  info1 ->SetTextSize(0.048);
  sprintf(legname,"#splitline{#chi^{2} = %4.2f}{dof = %i}",f1 -> GetChisquare(),f1 -> GetNDF());
  info1->DrawLatex(0.22,0.82,legname);


  AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta-1],etaBins[eta]);
  info1->DrawLatex(0.6,0.7,AuxString);

  pdfFile = (TString) "plots/Ratio_Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_data_comparison_" + method + (TString) ".pdf";
  c11 -> SaveAs(pdfFile);
    
  return 0;
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotDataMCComparisonFINAL(){

  setTDRStyle(true);
  gROOT->ForceStyle();
  
  const TString JetType = "PFCHS";
  TString Method;
  
  const int nEta =nEtaBins;

  TString etaString, filename;   
  

  TString rootFiles, AuxString;
  TGraphErrors* ratioEtaBinned[3] = {0};
 

  for(int meth = 0; meth<3; meth ++){

    double *ratioEtaBinnedX  = new double[nEta];
    double *ratioEtaBinnedY  = new double[nEta];
    double *ratioEtaBinnedEX = new double[nEta];
    double *ratioEtaBinnedEY = new double[nEta];

    if(meth == 0)      Method  = "RMS95";
    else if(meth == 1) Method  = "RMS99";
    else if(meth == 2) Method  = "gaus";


    cout<<endl<<endl<<"METHOD :   "<<Method<<endl;
  
  for(int eta = 0; eta < nEtaBins; eta++){
    
    cout<< endl<<endl<<endl<<eta+1<<". eta Bin!!"<<endl;
    
    // Read the MC and data results 
    rootFiles = (TString) "../plots_2012/PF_L1CHS/data/root_files/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + JetType + (TString) "_data_" + Method + (TString) ".root";
    TGraphErrors* JERData = readTGraphErrors(rootFiles,"Graph;1","Graph;1");
    rootFiles = (TString) "../plots_2012/PF_L1CHS/mc/root_files/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + JetType + (TString) "_mc_" + Method + (TString) ".root";
    TGraphErrors* JERMC = readTGraphErrors(rootFiles,"Graph","Graph");
    
    if(eta+1 == 1) etaString = Form("JER for |#eta| < %4.1f",etaBins[eta+1]);
    else           etaString = Form("JER for %4.1f <|#eta|< %4.1f",etaBins[eta+1],etaBins[eta+2]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // 1.) Calculate the ratio w/o systematic Uncertainties  

    int nData    = JERData->GetN();
    cout<<"nData = "<<nData<<endl;

    double *dataX  = JERData->GetX();
    double *dataY  = JERData->GetY();
    double *dataEX = JERData->GetEX();
    double *dataEY = JERData->GetEY();
    
    double *mcX  = JERMC->GetX();
    double *mcY  = JERMC->GetY();
    double *mcEX = JERMC->GetEX();
    double *mcEY = JERMC->GetEY();
    
    double *ratioX  = new double[nData];
    double *ratioY  = new double[nData];
    double *ratioEX = new double[nData];
    double *ratioEY = new double[nData];

 
    int idxMC = 0;
    for(int i=0; i<nData; i++){


      if(TMath::Abs(dataX[i]/mcX[idxMC] - 1.) > 0.05){
	i -= 1;
	idxMC += 1;
	continue;
      }
     
      ratioX[i]  = 1./2.*(dataX[i] + mcX[idxMC]);
      ratioY[i]  = dataY[i]/mcY[idxMC];
      ratioEX[i] = 1./2.*TMath::Sqrt(TMath::Power(dataEX[i],2)+TMath::Power(mcEX[idxMC],2));
      ratioEY[i] = TMath::Sqrt(TMath::Power((1./mcY[idxMC]),2)*TMath::Power(dataEY[i],2)+TMath::Power((dataY[i]/(TMath::Power(mcY[idxMC],2))),2)*TMath::Power(mcEY[idxMC],2));
    }
    
    TGraphErrors *Ratio = new TGraphErrors(nData,ratioX,ratioY,ratioEX,ratioEY);

    if(eta+1 == 1 ) AuxString = Form("Ratio between Data and MC for |#eta| < %4.1f",etaBins[eta+1]);
    else            AuxString = Form("Ratio between Data and MC for %4.1f <|#eta|<%4.1f",etaBins[eta+1],etaBins[eta+2]);
 
    Ratio -> SetTitle(AuxString); 
    Ratio -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
    Ratio -> GetXaxis() -> SetTitleOffset(1.1); 
    Ratio -> GetYaxis() -> SetTitle("Ratio of JER (DATA/MC)");
    Ratio -> GetYaxis() -> SetTitleOffset(1.2);   
    Ratio -> GetXaxis() -> SetLimits(0,600);
    TF1* f1 = new TF1("name","pol0",0,600);   
    Ratio   -> Fit("name","QR");
   
   
    double fitPar      = f1 -> GetParameter(0);
    double fitParStatE = f1 -> GetParError(0);

    // Make another TGraph with the fitted error bands   

    double x[100000] = {0};
    double ex[100000] = {0};
    double y[100000] = {0};
    double ey[100000]= {0};

    for(int i=0; i<100000; i++){

      x[i] = 0. + i*1./999.*600.;
      ex[i] = 0.;
      y[i] = fitPar;
      ey[i] = fitParStatE;
    }


    TGraphErrors *band = new TGraphErrors(100000,x,y,ex,ey);

    cout<<endl<<endl<<"fitPar +- fitParStatE = "<<fixed<<setprecision(3)<<fitPar<<"+/-"<<fitParStatE<<endl<<endl;

   
    TCanvas *c11 = new TCanvas("c11",AuxString,200,10,500,500);
    c11 -> cd();
    Ratio -> SetMinimum(0.5);
    Ratio -> SetMaximum(2.0);
  
    Ratio  -> DrawClone("AP"); 
    band->SetLineColor(2);
    band->SetFillColor(2);
    band->SetFillStyle(3001);
    band->SetMarkerColor(2);
    band->SetMarkerSize(0.00001);
    band ->SetMinimum(0.5);
    band ->SetMaximum(0.5);
    band ->Draw("e3same");
    f1->SetLineColor(2);
    f1 ->SetLineWidth(1);
    f1 ->Draw("same");
    Ratio  -> Draw("Psame"); 
  
    TLatex*  info   = new TLatex();
    info -> SetTextFont(132);
    info -> SetNDC();
    info -> SetTextSize(0.05); 
    info -> DrawLatex(0.20,0.78,Form("#chi^{2}/ndof = %4.2f / %i",f1 -> GetChisquare(),f1 -> GetNDF()));
    info -> DrawLatex(0.60,0.78,Form("f = %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0)));
    info -> DrawLatex(0.60,0.85,Form("#bf{%4.1f < |#eta^{Jet}| < %4.1f}",etaBins[eta],etaBins[eta+1]));

    filename = (TString) "plots/Ratio_Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + JetType + (TString) "_data_comparison_" + Method + (TString) ".pdf";
    c11 -> SaveAs(filename);
    
    ratioEtaBinnedX[eta]  = (etaBins[eta+1] + etaBins[eta])/2.; 
    ratioEtaBinnedY[eta]  = f1 -> GetParameter(0);
    ratioEtaBinnedEX[eta] = 0;
    ratioEtaBinnedEY[eta] = f1->GetParError(0);
 
    Ratio -> SetMarkerColor(1);
    Ratio -> SetLineColor(1);
    Ratio -> SetMarkerStyle(20);
    Ratio -> GetFunction("name")->SetLineColor(1);
    Ratio -> Draw("AP");
  
    
    TLegend* legend = new TLegend(0.4,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(Ratio,"Central Value","l");
    legend -> Draw("same");
    if(meth == 2) cout<<"ratioEtaBinnedX = "<<ratioEtaBinnedX[eta]<<endl;
    if(meth == 2) cout<<"ratioEtaBinnedY = "<<ratioEtaBinnedY[eta]<<endl;
   
  }
  

  ratioEtaBinned[meth] = new TGraphErrors(nEta,ratioEtaBinnedX,ratioEtaBinnedY,ratioEtaBinnedEX,ratioEtaBinnedEY);

  
  filename = (TString) "plots/RatioEtaBinned_" + JetType + (TString) "_" + Method + (TString) ".root";
  TFile *f = new TFile(filename,"RECREATE");
  f -> WriteTObject(ratioEtaBinned[meth],"Graph");
  f->Close();
  delete f;

  delete ratioEtaBinnedX;
  delete ratioEtaBinnedY;
  delete ratioEtaBinnedEX;
  delete ratioEtaBinnedEY;

  }


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Comparison to 2011 Data with all three methods
  cout<<endl; 
  
  TCanvas *cFinal = new TCanvas("cFinal","cFinal",200,10,500,500);
  cFinal -> cd();  


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

  ratioEtaBinned[2] -> GetXaxis() -> SetTitle("|#eta|");
  ratioEtaBinned[2] -> GetXaxis() -> SetLimits(0., 2.3);
  ratioEtaBinned[2] -> GetXaxis() -> SetNdivisions(505, "X");
  ratioEtaBinned[2] -> GetYaxis() -> SetTitle("Data/MC ratio (const fit)");
  ratioEtaBinned[2] -> GetYaxis() -> SetRangeUser(0.8, 1.5);
  ratioEtaBinned[2] -> SetMarkerStyle(21); 
  ratioEtaBinned[2] -> SetMarkerSize(1.2);
  ratioEtaBinned[2] -> SetLineColor(kRed);
  ratioEtaBinned[2] -> SetMarkerColor(kRed);
  ratioEtaBinned[2] -> Draw("Ae1p");
  
  Res_2011->SetMarkerStyle(20);
  Res_2011->SetMarkerSize(1.4);
  Res_2011->SetFillColor(kGray);
  Res_2011->SetFillStyle(3001);
  Res_2011->SetLineColor(kGray);
  Res_2011->DrawClone("e3psame");


  ratioEtaBinned[1] -> SetMarkerStyle(23);
  ratioEtaBinned[1] -> SetMarkerSize(1.2);
  ratioEtaBinned[1] -> SetLineColor(kBlue-4);
  ratioEtaBinned[1] -> SetMarkerColor(kBlue-4);
  ratioEtaBinned[1] -> Draw("e1psame");

  ratioEtaBinned[0]->SetMarkerStyle(24);
  ratioEtaBinned[0]->SetMarkerSize(1.2);
  ratioEtaBinned[0]->SetLineColor(3);
  ratioEtaBinned[0]->SetMarkerColor(3);
  ratioEtaBinned[0]->Draw("e1psame");

  ratioEtaBinned[2] -> Draw("e1psame");

  Res_2011->SetPointError(0, 0., 0., 0., 0.);
  Res_2011->SetPointError(1, 0., 0., 0., 0.);
  Res_2011->SetPointError(2, 0., 0., 0., 0.);
  Res_2011->SetPointError(3, 0., 0., 0., 0.);
  Res_2011->SetPointError(4, 0., 0., 0., 0.);

  Res_2011->Draw("psame");

  cmsPrel();

  TLegend *leg = new TLegend(0.20, 0.70, 0.40, 0.90);
  leg->SetBorderSize(0);
  // leg->SetBorderMode(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  
  leg->AddEntry(Res_2011,"2011 (total error)", "pl");
  leg->AddEntry(ratioEtaBinned[2],"Gaus Fit 2 #sigma", "pl");
  leg->AddEntry(ratioEtaBinned[1],"RMS 99%", "pl");
  leg->AddEntry(ratioEtaBinned[0],"RMS 95%", "pl");
   
  leg->Draw("same");
  
  cFinal->Print("plots/resultsComparison.png");   
  cFinal->Print("plots/resultsComparison.pdf");
   

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
      cout<<"in"<<endl;
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



/*
  namespace util
  {
  class FileOps
  {
  public:
  static TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName = "");
  static TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName = "");
  static TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName = "");
  static TH2* readTH2(const TString &fileName, const TString &hName, const TString &newHName = "", bool useCurrentStyle = true);
  static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  static util::HistVec readHistVec(const TString &fileName, const std::vector<TString> &histNames, const TString &newHistNameSuffix = "", bool useCurrentStyle = true);
  static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  static std::vector<TF1*> readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName = "");
  static THStack* readTHStack(const TString &fileName, const TString &stackName, const TString &newStackName = "");
  static TChain* createTChain(const TString &fileName, const TString &treeName = "", unsigned int verbosity = 1);
  };


  // -------------------------------------------------------------------------------------
  TH2* FileOps::readTH2(const TString &fileName, const TString &hName, const TString &newHName, bool useCurrentStyle) {
  TFile file(fileName,"READ");
  TH2* h = 0;
  file.GetObject(hName,h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  if( newHName.Length() ) h->SetName(newHName);
  } else {
  std::cerr << "ERROR in FileOps::readTH2: No TH2 with name '" << hName << "' in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return h;
  }


  //! Read TH1 histogram from file
  // -------------------------------------------------------------------------------------
  TH1* FileOps::readTH1(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
  TFile file(fileName,"READ");
  TH1 *h = 0;
  file.GetObject(histName,h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  if( newHistName.Length() ) h->SetName(newHistName);
  } else {
  std::cerr << "ERROR in FileOps::readTH1: Histogram with name '" << histName << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return h;
  }


  //! Read THStack from file
  // -------------------------------------------------------------------------------------
  THStack* FileOps::readTHStack(const TString &fileName, const TString &stackName, const TString &newStackName) {
  TFile file(fileName,"READ");
  THStack *s = 0;
  file.GetObject(stackName,s);
  if( s ) {
  //      s->SetDirectory(0);
  if( newStackName.Length() ) s->SetName(newStackName);
  } else {
  std::cerr << "ERROR in FileOps::readTHStack: THStack with name '" << stackName << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return s;
  }



  //! Read TGraph from file
  // -------------------------------------------------------------------------------------
  TGraph* FileOps::readTGraph(const TString &fileName, const TString &gName, const TString &newGName) {
  TFile file(fileName,"READ");
  TGraph *g = 0;
  file.GetObject(gName,g);
  if( g ) {
  if( newGName.Length() ) g->SetName(newGName);
  } else {
  std::cerr << "ERROR in FileOps::readTGraph: TGraph with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  file.Close();
    
  return g;
  }
  
  
  
  
  
  //! Read TH1 histograms from different files
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
  util::HistVec v(fileNames.size());
  for(unsigned int i = 0; i < fileNames.size(); ++i) {
  v[i] = readTH1(fileNames[i],histName,newHistName+util::toTString(i),useCurrentStyle);
  }
  std::cout << "Done\n";
    
  return v;
  }
  

  //! Read TH1 histograms from one file
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readHistVec(const TString &fileName, const std::vector<TString> &histNames, const TString &newHistNameSuffix, bool useCurrentStyle) {
  util::HistVec v;
  TFile file(fileName,"READ");
  for(std::vector<TString>::const_iterator it = histNames.begin();
  it != histNames.end(); ++it) {
  TH1 *h = 0;
  file.GetObject(*it,h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  h->SetName((*it)+newHistNameSuffix);
  v.push_back(h);
  } else {
  std::cerr << "ERROR in FileOps::readHistVec: Histogram with name '" << *it << "' does not exist in file '" << fileName << "'\n.";
  file.Close();
  exit(-1);
  }
  }
  file.Close();
    
  return v;
  }


  
  //! Read TH1 histograms from one file
  // -------------------------------------------------------------------------------------
  util::HistVec FileOps::readHistVec(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
  util::HistVec v;
  TFile file(fileName,"READ");
  bool binExists = true;
  unsigned int bin = 0;
  while( binExists ) {
  TH1 *h = 0;
  file.GetObject(histName+util::toTString(bin),h);
  if( h ) {
  h->SetDirectory(0);
  if( useCurrentStyle ) h->UseCurrentStyle();
  if( newHistName.Length() ) h->SetName(newHistName+util::toTString(bin));
  v.push_back(h);
  ++bin;
  } else {
  binExists = false;
  }
  }
  file.Close();
    
  if( v.size() == 0 ) std::cerr << "WARNING in util::FileOps::readHistVec(): No histogram read!\n";
    
  return v;
  }


  // -------------------------------------------------------------------------------------
  TF1* FileOps::readTF1(const TString &fileName, const TString &fName, const TString &newFName) {
  TFile file(fileName,"READ");
  TF1 *f = 0;
  file.GetObject(fName,f);
  if( f ) {
  //f->SetDirectory(0);
  if( newFName.Length() ) f->SetName(newFName);
  } else {
  std::cerr << "ERROR in FileOps::readTF1: TF1 with name '" << fName << "' does not exist in file '" << fileName << "'\n.";
  }
  file.Close();
    
  return f;
  }


  // -------------------------------------------------------------------------------------
  std::vector<TF1*> FileOps::readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName) {
  std::vector<TF1*> v;
  TFile file(fileName,"READ");
  bool binExists = true;
  unsigned int bin = 0;
  while( binExists ) {
  TF1* f = 0;
  file.GetObject(fName+util::toTString(bin),f);
  if( f ) {
  //f->SetDirectory(0);
  if( newFName.Length() ) f->SetName(newFName+util::toTString(bin));
  v.push_back(f);
  ++bin;
  } else {
  binExists = false;
  }
  }
  file.Close();
    
  if( v.size() == 0 ) std::cerr << "WARNING in util::FileOps::readTF1Vec(): No TF1 read!\n";
    
  return v;
  }



  //! Create TChain from input root files. The root
  //! files are expected to contain a TTree "DiJetTree".
  //! There are two possible input options:
  //!
  //! 1) 'fileName' specifies a single root file; it ends
  //!    with '.root';
  //! 2) 'fileName' contains a list of root file names.
  // --------------------------------------------------
  TChain* FileOps::createTChain(const TString &fileName, const TString &treeName, unsigned int verbosity) {
  if( verbosity >= 1 ) std::cout << "Creating TChain" << std::endl;

  TString tree = (treeName=="" ? "DiJetTree" : treeName);
  TChain* chain = new TChain(tree); 
    
  // Option 1: single root file
  if( fileName.EndsWith(".root") ) {
  if( verbosity >= 1 ) std::cout << "  Adding '" << tree << "' from file '" << fileName << "'" << std::endl;
  chain->Add(fileName);
  }
  // Option 2: list of root files
  else {
  if( verbosity >= 1 ) std::cout << "  Opening files from list '" << fileName << "'" << std::endl;
  std::ifstream filelist;
  filelist.open(fileName);
  int nOpenedFiles = 0;
  if( filelist.is_open() ) {
  TString name = "";
  while( !filelist.eof() ) {
  filelist >> name;
  if( filelist.eof() ) break;
  if( verbosity >= 1 ) std::cout << "  Adding '" << tree << "' from file '" << name << "'" << std::endl;
  chain->Add(name);
  nOpenedFiles++;
  }
  } else {
  std::cerr << "ERROR opening file '" << fileName << "'\n";
  exit(1);
  }
  filelist.close();
  }
    
  return chain;
  }
  }
*/







