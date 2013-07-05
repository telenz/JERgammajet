//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------- Script to create the MC Pileup distribution -------------------------------------------- 
//------------------------------------------------------------------------------------------------------------------------------------
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>

const int numTrigger = 8;

TH1D* readTH1(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle); 
TH1D* createTH1(TH1D *histo,const TString &HistoName,const TString &HistoTitle,int nBin, double bound_low, double bound_up,const TString &xTitle);
void plotTH1(TH1* histo,const TString &filename);

void calcPUMCdist(int nBins){
  
  // ----- Create Histogram for PU distribution -----
  //TH1D *hPUMCDist = new TH1D("hPUMCDist","PU Distribution in MC",600,0,60);
  TH1D *hPUMCDist = new TH1D("hPUMCDist","PU Distribution in MC",nBins,0,60);
  hPUMCDist -> Sumw2();
  hPUMCDist -> SetXTitle("# PU");
  
  // ----- Add Chain to MC Sample
  TString DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/pythia_flat_535_START53_V22/MCPhoton2012/";
  DataFilename += "ak5PFCHS.root";  
  TChain* chain = new TChain("GammaJetTree");
  chain -> Add(DataFilename); 

  // Plus QCD Sample
  //DataFilename = "/scratch/hh/dust/naf/cms/user/telenz/mc/qcd_pythia_flat_535_START53_V20/";
  //DataFilename += "ak5PFCHS*.root"; 
  //chain -> Add(DataFilename); 
  
  // ----- Set branch address to interesting variables
  float weight;
  float PUMCNumTruth; 
  int vtxN; 
  chain->SetBranchAddress("Weight",&weight);
  chain->SetBranchAddress("PUMCNumTruth",&PUMCNumTruth);
  chain->SetBranchAddress("VtxN",&vtxN);
  
  //  ----- Go through alle events and fill the MC Pileup-distribution with the true Pileup number
  cout<<"There are "<<chain->GetEntries()<<" entries in the Sample!!"<<endl<<endl;
  for(int n = 0; n < chain->GetEntries(); n++) {
    //for(int n = 0; n < 10000; n++) {
    chain->GetEntry(n);
    if( (n+1)%1000000 == 0 ) std::cout << "Event " << (n+1) << std::endl;
    hPUMCDist -> Fill(PUMCNumTruth,weight);
  }
  
  // ----- Normalize the histogram to one
  //hPUMCDist -> Scale(1.0/hPUMCDist->Integral());     
  //cout<<"Scale "
  
  // ----- Save and plot the PU distribution histogram
  TCanvas* canvas = new TCanvas("PUMC","PUMCDist",0,0,500,500);
  canvas -> SetLeftMargin(0.17);  
  hPUMCDist ->GetYaxis()->SetTitleOffset(1.5); 
  canvas -> cd();
  hPUMCDist -> Draw(); 
  canvas -> SaveAs(Form("cmssw5_3_5_Tag04_02_02LumiCalc2FineBinning_MBXS73000/weights_%iBins/PUMCdist.pdf",nBins));
  TFile *f = new TFile(Form("cmssw5_3_5_Tag04_02_02LumiCalc2FineBinning_MBXS73000/weights_%iBins/PUMCdist.pdf",nBins),"RECREATE");
  f -> WriteTObject(hPUMCDist);
  f ->Close();
  delete f;
  

  // ----- Read Data PU Distribution and calculate weight histograms for each trigger
  
  // ----- Read Pileup-distribution from Data (seperately for different Triggers)
  TH1D* hPUDataDist[numTrigger] = {0};
  TH1D* hPUWeight[numTrigger]   = {0};
  TString histoName;
  TString title;
  for(int i=0; i<numTrigger; i++){
    title.Form("PU-Weights for %i. pT Bin",i);
    histoName.Form("PUWeight%i",i);
    hPUWeight[i] = createTH1(hPUWeight[i],histoName,title,60,0,60,"# PU");
  }
  
  TString PUVersion = "cmssw5_3_5_Tag04_02_02LumiCalc2FineBinning_MBXS73000/";
  
  hPUDataDist[0] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon20_CaloIdVL_IsoL.root","pileup;1","pileup0", 0);
  hPUDataDist[1] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon30_CaloIdVL_IsoL.root","pileup;1","pileup1", 0);
  hPUDataDist[2] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon50_CaloIdVL_IsoL.root","pileup;1","pileup2", 0);
  hPUDataDist[3] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon75_CaloIdVL_IsoL.root","pileup;1","pileup3", 0);
  hPUDataDist[4] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon90_CaloIdVL_IsoL.root","pileup;1","pileup4", 0);
  hPUDataDist[5] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon135.root","pileup;1","pileup5", 0);
  hPUDataDist[6] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon150.root","pileup;1","pileup6", 0);
  hPUDataDist[7] =  readTH1(PUVersion + "MyDataPileupHistogramHLT_Photon160.root","pileup;1","pileup7", 0);

  if(nBins == 60) for(int i = 0; i<numTrigger; i++ ) hPUDataDist[i] -> Rebin(10);
  
  // ----- Normalize PUData Distribution to 1.0
  for(int i = 0; i<numTrigger; i++ ) hPUDataDist[i] -> Scale(hPUMCDist->Integral()/hPUDataDist[i]->Integral());
  
  // ----- Calculate weights and normalize weight histogram to 1
  for(int i=0;i<numTrigger;i++) {
    cout<<hPUDataDist[i]->GetNbinsX()<<endl;
    cout<<hPUMCDist->GetNbinsX()<<endl;
    hPUWeight[i] -> Divide(hPUDataDist[i],hPUMCDist);   
    //hPUWeight[i] -> Scale(1.0/hPUWeight[i]->Integral());          
  }
    
  // ----- Plot the weights for all different pTGamma bins  
  for(int i=0; i<numTrigger; i++){
    title.Form("cmssw5_3_5_Tag04_02_02LumiCalc2FineBinning_MBXS73000/weights_%iBins/PUWeights_%i",nBins,i);
    plotTH1(hPUWeight[i],title);
  }
   
}

//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------- Script to create the MC Pileup distribution -------------------------------------------- 
//------------------------------------------------------------------------------------------------------------------------------------


TH1D* readTH1(const TString &fileName, const TString &histName, const TString &newHistName, bool useCurrentStyle) {
 
  TFile f(fileName,"READ");
  TH1D *h = 0;
  f.GetObject(histName,h);
  if( h ) {
    h->SetDirectory(0);
    if( useCurrentStyle ) h->UseCurrentStyle();
    if( newHistName.Length() ) h->SetName(newHistName);
  } else {
    std::cerr << "ERROR in FileOps::readTH1: Histogram with name '" << histName << "' does not exist in file '" << fileName << "'\n.";
    f.Close();
    exit(-1);
  }
  f.Close();
  
  return h;
}


void plotTH1(TH1* histo,const TString &filename){

  TString tot_filename;
  TFile* f;

  TCanvas* canvas = new TCanvas(filename,filename,0,0,500,500);
  
  canvas -> SetLeftMargin(0.17);
 
  histo ->GetYaxis()->SetTitleOffset(1.5); 
  canvas -> cd();
  histo -> Draw(); 
  
  
  tot_filename = filename + ".pdf";
  canvas -> SaveAs(tot_filename);
  
  tot_filename = filename + ".root";
  f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(histo);
  f ->Close();
  delete f;
}


TH1D* createTH1(TH1D *histo,const TString &HistoName,const TString &HistoTitle,int nBin, double bound_low, double bound_up,const TString &xTitle){
  
  histo = new TH1D(HistoName,HistoTitle,nBin,bound_low,bound_up);
  histo -> Sumw2();
  histo -> SetXTitle(xTitle);
    
  return histo;
}
//------------------------------------------------------------------------------------------------------------------------------------
