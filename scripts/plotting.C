// This file is supposed to plot TGraphs also more than 1 in one plot

#ifndef FileOps_h
#define FileOps_h

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
#include "TText.h"
#include "TROOT.h"


//#include "utils.h"
//#include "HistOps.h"

TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName);
//TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName);
TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);
TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName);

const int pt_int = 12;
const int eta_int = 4;

int plotJetResponseAlpha(int variable){


  gStyle -> SetTitleSize(0.05,"X"); 
  gStyle -> SetTitleSize(0.05,"Y"); 
  gStyle -> SetLabelSize(0.04,"X");
  gStyle -> SetLabelSize(0.04,"Y");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.3,"X");
  gStyle -> SetTitleOffset(1.3,"Y");
  //gStyle -> SetOptLogy();
  gStyle -> SetHistLineWidth(1);
  gROOT->ForceStyle(); 
  gStyle->SetOptStat("");


  gDirectory->Delete(); 
  gROOT->GetListOfCanvases()->Delete();

  TString sourceResponseAlpha, sourceIntrinsicAlpha, sourceImbalanceAlpha, sourceTotalAlpha, pdfFile, etaRegion, ptRegion, title, EtaPtRegion; 
  double ptBorders[pt_int+1]   = {0}; 
  double etaBorders[eta_int+1] = {0};

  etaBorders[0] = 0.0;
  etaBorders[1] = 0.5;
  etaBorders[2] = 1.1; 
  etaBorders[3] = 1.7; 
  etaBorders[4] = 2.3; 

  ptBorders[0]  = 22;
  ptBorders[1]  = 33;
  ptBorders[2]  = 55;
  ptBorders[3]  = 82.5;
  ptBorders[4]  = 99;
  ptBorders[5]  = 148.5;
  ptBorders[6]  = 165;
  ptBorders[7]  = 176;
  ptBorders[8]  = 200;
  ptBorders[9]  = 250;
  ptBorders[10] = 300;
  ptBorders[11] = 400;

  
  TLatex*  info;
  TLegend *legend;
  TGraphErrors* Add;
  TF1 *fScaleAlpha;
  TF1 *total;
  for(int j=0; j<eta_int; j++){
    for(int i=0;i<pt_int;i++){       

      if(variable ==1){
	//sourceResponseAlpha.Form("plots_2012/mc/PF_L1FastJet/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_PF_mc.root",j+1,i+1);
	//sourceIntrinsicAlpha.Form("plots_2012/mc/PF_L1FastJet/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_intrinsic_PF_mc.root",j+1,i+1);
	//sourceImbalanceAlpha.Form("plots_2012/mc/PF_L1FastJet/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_imbalance_PF_mc.root",j+1,i+1);
	//sourceTotalAlpha.Form("plots_2012/mc/PF_L1FastJet/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_total_PF_mc.root",j+1,i+1);
	//pdfFile.Form("plots_2012/PF_L1FastJet/mc/pdfs/JES_for_%i_eta_bin_%i_pTGamma_bin_ALL_PF_mc.pdf",j+1,i+1); 
	sourceResponseAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_PF_mc.root",j+1,i+1);
	sourceIntrinsicAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_intrinsic_PF_mc.root",j+1,i+1);
	sourceImbalanceAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_imbalance_PF_mc.root",j+1,i+1);
	sourceTotalAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_scale_for_%i_eta_bin_%i_pTGamma_bin_total_PF_mc.root",j+1,i+1);
	pdfFile.Form("plots_2012/PF_L1FastJet/mc/pdfs/JES_for_%i_eta_bin_%i_pTGamma_bin_ALL_PF_mc.pdf",j+1,i+1); 
	title = "Jet Energy Scale";
  
      }
      else if(variable == 2){
	//sourceResponseAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_PF_mc.root",j+1,i+1);
	//sourceIntrinsicAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_intrinsic_PF_mc.root",j+1,i+1);
	//sourceImbalanceAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_imbalance_PF_mc.root",j+1,i+1);
	//sourceTotalAlpha.Form("plots_2012/PF_L1FastJet/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_total_PF_mc.root",j+1,i+1);

	sourceResponseAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_PFCHS_mc_RMS99.root",j+1,i+1);
	sourceIntrinsicAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_intrinsic_PFCHS_mc_RMS99.root",j+1,i+1);
	sourceImbalanceAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_imbalance_PFCHS_mc_RMS99.root",j+1,i+1);
	sourceTotalAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_total_PFCHS_mc_RMS99.root",j+1,i+1);
	//pdfFile.Form("plots_2012/PF_L1FastJet/mc/JER_for_%i_eta_bin_%i_pTGamma_bin_ALL_PF_mc.pdf",j+1,i+1);
	pdfFile.Form("plots/JER_for_%i_eta_bin_%i_pTGamma_bin_all_contributions_PFCHS_RMS99_mc.pdf",j+1,i+1);
	title = "Jet Energy Resolution";
      }
  
 
  
      

      TCanvas *c = new TCanvas(EtaPtRegion,title,200,10,450,450);
      c -> cd();
      
      legend  = new TLegend(0.45,0.7,0.9,0.9);
      legend -> SetFillColor(0);
      legend -> SetTextFont(132);
      legend -> SetTextSize(0.042);
  
      Add = readTGraphErrors(sourceResponseAlpha,"Graph;1","Graph;1");  
      if(Add == 0) continue;
      
      if(variable == 1) fScaleAlpha = Add->GetFunction("fScaleAlpha");
      else              fScaleAlpha = Add->GetFunction("fResolutionAlpha");
      Add -> SetMarkerColor(1);
      Add -> SetLineColor(1);
      fScaleAlpha -> SetLineColor(1);
      Add -> GetXaxis() -> SetTitle("#alpha [%]");   
      Add -> GetXaxis() -> SetTitle("p_{T}^{2nd jet}/p_{T}^{#gamma} #upoint 100");    
      Add -> GetYaxis() -> SetTitle("Resolution"); 
      Add -> SetTitle("");
            
      
      legend -> AddEntry(Add,"#gamma + Jet (pseudo data)","l");
      Add -> Draw("AP");

      
      if(variable ==1){
	Add -> GetYaxis() -> SetTitle("JES");   
	Add -> SetMinimum(0.8);
	Add -> SetMaximum(1.1);  
      }
      else if(variable ==2){
	//Add -> GetYaxis() -> SetTitle("JER");
	Add -> SetMinimum(0.0);
	Add -> SetMaximum(0.4);   
      }
  
      Add -> GetXaxis()->SetLimits(0,20);
      
  
      // Intrinsic
      Add = readTGraphErrors(sourceIntrinsicAlpha,"Graph;1","Graph;1");  
      if(variable == 1) fScaleAlpha = Add->GetFunction("fScaleAlpha");
      else              fScaleAlpha = Add->GetFunction("fResolutionAlpha");
  
      Add -> SetMarkerColor(46);
      Add -> SetLineColor(46);
      fScaleAlpha -> SetLineColor(46);

      legend -> AddEntry(Add,"Intrinsic","l");
      Add -> Draw("Psame");

      Add = readTGraphErrors(sourceImbalanceAlpha,"Graph;1","Graph;1");  
      if(variable == 1) fScaleAlpha = Add->GetFunction("fScaleAlpha");
      else              fScaleAlpha = Add->GetFunction("fResolutionAlpha");
  
      Add -> SetMarkerColor(4);
      Add -> SetLineColor(4);
      fScaleAlpha -> SetLineColor(4);
      legend -> AddEntry(Add,"Imbalance","l");
      Add -> Draw("Psame");
      
      if(variable == 1) total = readTF1(sourceTotalAlpha,"totalScale;1", "totalScale;1");
      else              total = readTF1(sourceTotalAlpha,"totalResolution;1", "totalResolution;1");
      total -> SetLineColor(14);
      total -> SetLineWidth(3);
      legend -> AddEntry(total,"Total","l");
      total -> Draw("same");

      // Draw info boxes
      if (i==pt_int-1) ptRegion.Form(" %4.1f GeV < p_{T}^{#gamma} \n",ptBorders[i]);
      else       ptRegion.Form("%4.1f GeV < p_{T}^{#gamma} < %4.1f GeV \n" ,ptBorders[i], ptBorders[i+1]);
      
      if(j == 0) etaRegion.Form("|#eta^{jet}| < %4.1f", etaBorders[j+1]);
      else       etaRegion.Form("%4.1f < |#eta^{jet}| < %4.1f",etaBorders[j], etaBorders[j+1]);
      
      EtaPtRegion.Form("#splitline{%s}{%s}",ptRegion.Data(),etaRegion.Data());

      info   = new TLatex();
      info -> SetTextFont(132);
      info -> SetNDC();
      info -> SetTextSize(0.042);
      if(variable == 1){
	info->DrawLatex(0.22,0.68,  ptRegion);
	info->DrawLatex(0.22,0.56,  etaRegion);
	info->DrawLatex(0.55,0.56, "Anti-k_{T} 0.5 PFCHSJets");
	
      }
      else{
	info->DrawLatex(0.20,0.63,  ptRegion);
	info->DrawLatex(0.20,0.56,  etaRegion);
	info->DrawLatex(0.50,0.56, "Anti-k_{T} 0.5 PFCHSJets");
      }

      
      legend -> Draw("same");
      c -> SaveAs(pdfFile);
    }
  }

  
  
  
  

  // FINAL PLOTS
  
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
    //exit(-1);
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


// -------------------------------------------------------------------------------------
TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName) {
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

//! Read TGraphAsymmErrors from file
// -------------------------------------------------------------------------------------
/*  TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName) {
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
    }*/



/*
  namespace util
  {
  class FileOps
  {
  public:
  static TF1* readTF1(const TString &fileName, const TString &fName, const TString &newFName = "");
  static TGraph* readTGraph(const TString &fileName, const TString &gName, const TString &newGName = "");
  static TGraphAsymmErrors* readTGraphAsymmErrors(const TString &fileName, const TString &gName, const TString &newGName = "");
  //static TH2* readTH2(const TString &fileName, const TString &hName, const TString &newHName = "", bool useCurrentStyle = true);
  //static TH1* readTH1(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  //static util::HistVec readTH1(const std::vector<TString> &fileNames, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  //static util::HistVec readHistVec(const TString &fileName, const std::vector<TString> &histNames, const TString &newHistNameSuffix = "", bool useCurrentStyle = true);
  //static util::HistVec readHistVec(const TString &fileName, const TString &histName, const TString &newHistName = "", bool useCurrentStyle = true);
  //static std::vector<TF1*> readTF1Vec(const TString &fileName, const TString &fName, const TString &newFName = "");
  //static THStack* readTHStack(const TString &fileName, const TString &stackName, const TString &newStackName = "");
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
#endif





