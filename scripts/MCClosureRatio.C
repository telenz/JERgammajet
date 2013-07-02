// ------------------------------------------------------------------------------------------------
// -------  Script to compare Jet Energy Resolutions and Scales between Data and MC (T.L.)  ------- 
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
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
#include "../CODE/myDeclarations.h"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);



int plotDataMCComparison(int eta, int variable){


  const TString method = "RMS99";

  //Comparison and Ratio Plots
  TString etaString; 
  char rootFileMC[100];
  char rootFileData[100];
  char pdfFile[100];
    
  if(variable ==1){
    sprintf(rootFileMC,"Scale_for_%d_eta_bin_PF_data_2011.root",eta);
    sprintf(rootFileData,"Scale_for_%d_eta_bin_PF_data_2012.root",eta);
    sprintf(pdfFile,"Scale_for_%d_eta_bin_PF_data_comparison.pdf",eta);
    if(eta == 1) etaString = "Jet Energy Scale for |eta|<0.5";
    else if(eta == 2) etaString = "Jet Energy Scale for 0.5<|eta|<1.3";
  }
  else if(variable == 2){

    // Fast Jets
    //sprintf(rootFileMC,"../jetphoton/plots_2012/PF_L1FastJet/mc/root_files/Resolution_for_%d_eta_bin_PF_mc.root",etaBins);
    //sprintf(rootFileData,"../jetphoton/plots_2012/PF_L1FastJet/data/root_files/Resolution_for_%d_eta_bin_PF_data.root",etaBins);
    // CHS Jets

    if(method == "RMS95"){
    sprintf(rootFileMC,"../plots_2012/PF_L1CHS/mc/root_files_NOTSMEARED_NEW/Resolution_for_%d_eta_bin_PFCHS_mc_RMS95.root",eta);
    sprintf(rootFileData,"../plots_2012/PF_L1CHS/mc/root_files_SMEARED_NEW/Resolution_for_%d_eta_bin_PFCHS_mc_RMS95.root",eta);

    //sprintf(rootFileMC,"../plots_2012/PF_L1CHS/mc/root_files_woPixelSeedCut/Resolution_for_%d_eta_bin_PFCHS_mc_RMS95.root",eta);
    //sprintf(rootFileData,"../plots_2012/PF_L1CHS/mc/root_files_withPixelSeedCut/Resolution_for_%d_eta_bin_PFCHS_mc_RMS95.root",eta);
    }
    else if(method == "RMS99"){
    sprintf(rootFileMC,"../plots_2012/PF_L1CHS/mc/root_files_NOTSMEARED_NEW/Resolution_for_%d_eta_bin_PFCHS_mc_RMS99.root",eta);
    sprintf(rootFileData,"../plots_2012/PF_L1CHS/mc/root_files_SMEARED_NEW/Resolution_for_%d_eta_bin_PFCHS_mc_RMS99.root",eta);
    }
    else if(method == "gaus"){
    sprintf(rootFileMC,"../plots_2012/PF_L1CHS/mc/root_files_NOTSMEARED_NEW/Resolution_for_%d_eta_bin_PFCHS_mc_gaus.root",eta);
    sprintf(rootFileData,"../plots_2012/PF_L1CHS/mc/root_files_SMEARED_NEW/Resolution_for_%d_eta_bin_PFCHS_mc_gaus.root",eta);
    }
    //sprintf(rootFileMC,"../jetphoton/plots_2012/PF_L1CHS/mc/root_files_MC_2sigma/Resolution_for_%d_eta_bin_PFCHS_mc.root",etaBins);
    //sprintf(rootFileMC,"../jetphoton/plots_2012/PF_L1CHS/mc/root_files_MC_2sigma_200/Resolution_for_%d_eta_bin_PFCHS_mc.root",etaBins);
    //sprintf(rootFileMC,"../jetphoton/plots_2012/PF_L1CHS/mc/root_files_MC_RMS95_EvenFinerBinning/Resolution_for_%d_eta_bin_PFCHS_mc.root",etaBins);
    //sprintf(rootFileMC,"../jetphoton/plots_2012/PF_L1CHS/mc/root_files_MC_RMS95/Resolution_for_%d_eta_bin_PFCHS_mc.root",etaBins);
    //sprintf(rootFileData,"../jetphoton/plots_2012/PF_L1CHS/mc/root_files/Resolution_for_%d_eta_bin_PFCHS_mc.root",etaBins);

    //2011
    //sprintf(rootFileMC,"../jetphoton/plots_2011/PF_L1FastJet/mc/root_files/Resolution_for_%d_eta_bin_PF_mc.root",etaBins);
    //sprintf(rootFileData,"../jetphoton/plots_2011/PF_L1FastJet/data/root_files/Resolution_for_%d_eta_bin_PF_data.root",etaBins);
   
    
    

    if(eta == 1)      etaString = Form("Jet Energy Resolution for |#eta^{Jet}| < %4.1f",etaBins[1]);
    else if(eta == 2) etaString = Form("Jet Energy Resolution for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[1],etaBins[2]);
    else if(eta == 3) etaString = Form("Jet Energy Resolution for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[2],etaBins[3]);
    else if(eta == 4) etaString = Form("Jet Energy Resolution for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[3],etaBins[4]);
    else if(eta == 5) etaString = Form("Jet Energy Resolution for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[4],etaBins[5]);
  }
    
  TGraphErrors* JERMC = readTGraphErrors(rootFileMC,"Graph;1","Graph;1");  
  JERMC -> SetMarkerColor(3);
  JERMC -> SetMarkerStyle(20);
  JERMC -> SetLineColor(3);
  //TCanvas* c1 = new TCanvas("c1","c1",100,100,500,500);
  //c1->cd();
  //JERMC -> Draw("AP");
    
  TGraphErrors* JERData = readTGraphErrors(rootFileData,"Graph;1","Graph;1");
  JERData -> SetMarkerColor(4);
  JERData -> SetMarkerStyle(20);
  JERData -> SetLineColor(4);
  //TCanvas* c2 = new TCanvas("c2","c2",100,100,500,500);
  //c2->cd();
  //JERData -> Draw("AP");
  
  TCanvas *c = new TCanvas("c",etaString,200,10,800,800);
  c -> SetLeftMargin(0.17);
  c -> cd();
  
  TMultiGraph* mg = new TMultiGraph();
  mg -> SetTitle(etaString);

  char fitName[100];
  if(variable == 1){
    sprintf(fitName,"fScale") ;
    JERMC   -> GetFunction(fitName) -> SetLineColor(8);
    JERData -> GetFunction(fitName) -> SetLineColor(9);
  }

  if(variable == 2){
    sprintf(fitName,"fResolution") ;
    JERMC   -> GetFunction(fitName) -> SetLineColor(3);
    JERData -> GetFunction(fitName) -> SetLineColor(4);
  }

  TLegend *legend  = new TLegend(0.5,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.033);
  //legend -> AddEntry(JERData,"Data","l");
  //legend -> AddEntry(JERMC,"MC","l");
  legend -> AddEntry(JERData,"with PixelSeed Cut","l");
  legend -> AddEntry(JERMC,"w/o PixelSeed Cut","l");
  
  mg -> Add(JERMC);
  mg -> Add(JERData);
  mg -> Draw("AP");
  
  mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma}");
  mg -> GetYaxis() -> SetTitleOffset(2.0); 
  //mg -> GetYaxis()->SetLimits(-0.4,0.3); 
  //mg -> GetXaxis()->SetLimits(0.,600.);
  
  if(variable ==1){
    mg -> GetYaxis() -> SetTitle("JES");
    mg -> SetMinimum(0.8);
    mg -> SetMaximum(1.1);  
  }
  else if(variable ==2){
    mg -> GetYaxis() -> SetTitle("JER");
    mg -> SetMinimum(0.0);
    mg -> SetMaximum(0.2);   
  }
  legend -> Draw("same");

  TString pdfName = (TString)"plots/MCClosure_Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_" + method + (TString) ".pdf";
  c -> SaveAs(pdfName);

  // -------------------  Ratio Plot for Resolution  -------------------

  int numEntries = 9;


  double axMC[11] = {0};
  double ayMC[11] = {0};
  double axData[11] = {0};
  double ayData[11] = {0};
  double axcomp[11] = {0};
  double aycomp[11] = {0};

  double axMCError[11] = {0};
  double ayMCError[11] = {0};
  double axDataError[11] = {0};
  double ayDataError[11] = {0};
  double axcompError[11] = {0};
  double aycompError[11] = {0};

  double auxX[numEntries], auxXErr[numEntries], auxY[numEntries],auxYErr[numEntries];
  
  int offsetMCData = 0;

  for(int i=0; i<numEntries; i++){
    JERMC   -> GetPoint(i+offsetMCData,axMC[i],ayMC[i]);
    //JERMC   -> GetPoint(i,axMC[i],ayMC[i]);
    JERData -> GetPoint(i+0,axData[i],ayData[i]);

    axMCError[i] = JERMC   -> GetErrorX(i+offsetMCData);
    ayMCError[i] = JERMC   -> GetErrorY(i+offsetMCData);
    axDataError[i] = JERData -> GetErrorX(i+0);
    ayDataError[i] = JERData -> GetErrorY(i+0);

    axcomp[i] = 1./2.*(axData[i]+axMC[i]);
    aycomp[i] = ayData[i]/ayMC[i];

    axcompError[i] =  1./2.*TMath::Sqrt(TMath::Power(axDataError[i],2)+TMath::Power(axMCError[i],2));
    aycompError[i] = TMath::Sqrt(TMath::Power((1./ayMC[i]),2)*TMath::Power(ayDataError[i],2)+TMath::Power((ayData[i]/(TMath::Power(ayMC[i],2))),2)*TMath::Power(ayMCError[i],2));

  }

  TString ratioName;

  for(int i=0;i<numEntries;i++){

    auxX[i] = axcomp[i];
    auxY[i] = aycomp[i];
    auxXErr[i] = axcompError[i];
    auxYErr[i] = aycompError[i];


    cout<<"axMC["<<i<<"] = "<<axMC[i]<<endl;
    cout<<"ayMC["<<i<<"] = "<<ayMC[i]<<endl;
    cout<<"ayMCError["<<i<<"] = "<<ayMCError[i]<<endl;
    cout<<"axData["<<i<<"] = "<<axData[i]<<endl;
    cout<<"ayData["<<i<<"] = "<<ayData[i]<<endl;
    cout<<"ayDataError["<<i<<"] = "<<ayDataError[i]<<endl;
    cout<<"Ratio["<<i<<"] = "<<aycomp[i]<<endl<<endl;


  }

  //TGraphErrors *Ratio = new TGraphErrors(10,axcomp,aycomp,axcompError,aycompError);
  TGraphErrors *Ratio = new TGraphErrors(numEntries,auxX,auxY,auxXErr,auxYErr);
  Ratio -> GetXaxis()->SetLimits(0,600);
  //char test[100] = "Ratio between Data and MC for";
  char test[100] = "Ratio  with/without pixelSeed Cut for";
  /*
  if(eta == 1 )     ratioName = Form("Ratio between Data and MC for |#eta^{Jet}| < %4.1f",etaBins[1]);
  else if(eta == 2) ratioName = Form("Ratio between Data and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[1],etaBins[2]);
  else if(eta == 3) ratioName = Form("Ratio between Data and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[2],etaBins[3]);
  else if(eta == 4) ratioName = Form("Ratio between Data and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[3],etaBins[4]);
  else if(eta == 5) ratioName = Form("Ratio between Data and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[4],etaBins[5]);
  */
  
  if(eta == 1 )     ratioName = Form("%s |#eta^{Jet}| < %4.1f",test,etaBins[1]);
  else if(eta == 2) ratioName = Form("%s  %4.1f < |#eta^{Jet}| < %4.1f",test,etaBins[1],etaBins[2]);
  else if(eta == 3) ratioName = Form("%s  %4.1f < |#eta^{Jet}| < %4.1f",test,etaBins[2],etaBins[3]);
  else if(eta == 4) ratioName = Form("%s  %4.1f < |#eta^{Jet}| < %4.1f",test,etaBins[3],etaBins[4]);
  else if(eta == 5) ratioName = Form("%s  %4.1f < |#eta^{Jet}| < %4.1f",test,etaBins[4],etaBins[5]);

  /*
  if(eta == 1 )     ratioName = Form("Ratio between MC_{smeared} and MC for |#eta^{Jet}| < %4.1f",etaBins[1]);
  else if(eta == 2) ratioName = Form("Ratio between MC_{smeared} and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[1],etaBins[2]);
  else if(eta == 3) ratioName = Form("Ratio between MC_{smeared} and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[2],etaBins[3]);
  else if(eta == 4) ratioName = Form("Ratio between MC_{smeared} and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[3],etaBins[4]);
  else if(eta == 5) ratioName = Form("Ratio between MC_{smeared} and MC for  %4.1f < |#eta^{Jet}| < %4.1f",etaBins[4],etaBins[5]);
  */
  Ratio -> SetTitle(ratioName); 
  Ratio -> GetXaxis() -> SetTitle("Photon pT");
  Ratio -> GetXaxis() -> SetTitleOffset(1.1); 
  //Ratio -> GetYaxis() -> SetTitle("Ratio of JER (MC_{smeared}/MC)");
  Ratio -> GetYaxis() -> SetTitleOffset(1.2);   

  TF1* f1 = new TF1("name","pol0",0,600);   
  Ratio -> Fit("name","");
  cout<<"ChiSquare = "<<f1 -> GetChisquare()<<endl;
  legend  = 0;
  legend = new TLegend(0.55,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.04);
  char legname[100];
  sprintf(legname," %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0));
  legend -> SetHeader(legname);
  TCanvas *c11 = new TCanvas("c11",ratioName,200,10,800,800);
  c11 -> cd();
  Ratio -> SetMinimum(0.5);
  Ratio -> SetMaximum(2.0);
  

  // Draw info boxes
  Ratio  -> Draw("AP"); 
  legend -> Draw("same");
  TLatex*  info   = new TLatex();
  
  info->SetTextFont(132);
  info-> SetNDC();
  info->SetTextSize(0.040);
  sprintf(legname,"#splitline{#chi^{2} = %4.2f}{dof = %i}",f1 -> GetChisquare(),f1 -> GetNDF());
  info->DrawLatex(0.22,0.84,legname);


  //if(eta ==1 ) info->DrawLatex(0.58,0.75,"input ratio: 1.133");
  //else if(eta ==2 ) info->DrawLatex(0.58,0.75,"input ratio: 1.083");
  //else if(eta ==3 ) info->DrawLatex(0.58,0.75,"input ratio: 1.145");
  //else if(eta ==4 ) info->DrawLatex(0.58,0.75,"input ratio: 1.288");
  
  pdfName = (TString)"plots/MCClosure_Ratio_Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_" + method + (TString) ".pdf";
  c11 -> SaveAs(pdfName);
  
  //delete JERMC;
  //delete JERData;
  //delete f1;
  //delete Ratio;
  
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







