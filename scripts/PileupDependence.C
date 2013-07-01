
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
#include "TH2.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"

#include "../CODE/myDeclarations.h"

TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool isData);

int makeGraphData(){

 
  TGraphErrors * graph;
  double x[7],y[7],xErr[7],yErr[7];
  TH2D* histo2D;
  

  TFile *file = TFile::Open("PhotonPtNVtx_2d_PF_data.root");
  file->GetObject("Photon Pt against Number of Vertices",histo2D);

  TH1D * histoVtx = (TH1D*)histo2D->ProjectionY("histoVtx",0,1200);

  
  TH1D* histo1D = (TH1D*)histo2D->ProjectionX("histo1D",1,5);
  TH1D* histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histoVtxPart->SetBinContent(1,histoVtx->GetBinContent(1));
  histoVtxPart->SetBinContent(2,histoVtx->GetBinContent(2));
  histoVtxPart->SetBinContent(3,histoVtx->GetBinContent(3));
  histoVtxPart->SetBinContent(4,histoVtx->GetBinContent(4));
  histoVtxPart->SetBinContent(5,histoVtx->GetBinContent(5));
  y[0] = histo1D->GetMean();
  yErr[0] = histo1D->GetMeanError();  
  x[0] = histoVtxPart->GetMean();
  xErr[0] = histoVtxPart->GetMeanError();

  cout<<"x0 = "<<x[0]<<endl;
  cout<<"x0Err = "<<xErr[0]<<endl;

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",6,10);
  histoVtxPart->SetBinContent(6,histoVtx->GetBinContent(6));
  histoVtxPart->SetBinContent(7,histoVtx->GetBinContent(7));
  histoVtxPart->SetBinContent(8,histoVtx->GetBinContent(8));
  histoVtxPart->SetBinContent(9,histoVtx->GetBinContent(9));
  histoVtxPart->SetBinContent(10,histoVtx->GetBinContent(10));
  y[1] = histo1D->GetMean();
  yErr[1] = histo1D->GetMeanError();
  x[1] = histoVtxPart->GetMean();
  xErr[1] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",11,15);
  histoVtxPart->SetBinContent(11,histoVtx->GetBinContent(11));
  histoVtxPart->SetBinContent(12,histoVtx->GetBinContent(12));
  histoVtxPart->SetBinContent(13,histoVtx->GetBinContent(13));
  histoVtxPart->SetBinContent(14,histoVtx->GetBinContent(14));
  histoVtxPart->SetBinContent(15,histoVtx->GetBinContent(15));
  y[2] = histo1D->GetMean();
  yErr[2] = histo1D->GetMeanError();
  x[2] = histoVtxPart->GetMean();
  xErr[2] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",16,20);
  histoVtxPart->SetBinContent(16,histoVtx->GetBinContent(16));
  histoVtxPart->SetBinContent(17,histoVtx->GetBinContent(17));
  histoVtxPart->SetBinContent(18,histoVtx->GetBinContent(18));
  histoVtxPart->SetBinContent(19,histoVtx->GetBinContent(19));
  histoVtxPart->SetBinContent(20,histoVtx->GetBinContent(20));
  y[3] = histo1D->GetMean();
  yErr[3] = histo1D->GetMeanError();
  x[3] = histoVtxPart->GetMean();
  xErr[3] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",21,25);
  histoVtxPart->SetBinContent(21,histoVtx->GetBinContent(21));
  histoVtxPart->SetBinContent(22,histoVtx->GetBinContent(22));
  histoVtxPart->SetBinContent(23,histoVtx->GetBinContent(23));
  histoVtxPart->SetBinContent(24,histoVtx->GetBinContent(24));
  histoVtxPart->SetBinContent(25,histoVtx->GetBinContent(25));
  y[4] = histo1D->GetMean();
  yErr[4] = histo1D->GetMeanError();
  x[4] = histoVtxPart->GetMean();
  xErr[4] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",26,30);
  histoVtxPart->SetBinContent(26,histoVtx->GetBinContent(26));
  histoVtxPart->SetBinContent(27,histoVtx->GetBinContent(27));
  histoVtxPart->SetBinContent(28,histoVtx->GetBinContent(28));
  histoVtxPart->SetBinContent(29,histoVtx->GetBinContent(29));
  histoVtxPart->SetBinContent(30,histoVtx->GetBinContent(30));
  y[5] = histo1D->GetMean();
  yErr[5] = histo1D->GetMeanError();
  x[5] = histoVtxPart->GetMean();
  xErr[5] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",31,35);
  histoVtxPart->SetBinContent(31,histoVtx->GetBinContent(31));
  histoVtxPart->SetBinContent(32,histoVtx->GetBinContent(32));
  histoVtxPart->SetBinContent(33,histoVtx->GetBinContent(33));
  histoVtxPart->SetBinContent(34,histoVtx->GetBinContent(34));
  histoVtxPart->SetBinContent(35,histoVtx->GetBinContent(35));
  y[6] = histo1D->GetMean();
  yErr[6] = histo1D->GetMeanError();
  x[6] = histoVtxPart->GetMean();
  xErr[6] = histoVtxPart->GetMeanError();
 

  graph = new TGraphErrors(7,x,y,xErr,yErr);
  TF1* f1 = new TF1("pol1","pol1"); 
  graph -> Fit(f1,"Q","",5.0,35.0);

  char legEntry[100];


  cout<<"par0 = "<<f1->GetParameter(0)<<endl;
  cout<<"par1 = "<<f1->GetParameter(1)<<endl;
  sprintf(legEntry,"%4.2f + %4.2f (#pm %4.2f) * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParError(1));

  TLegend* leg_hist       =  0 ;
  leg_hist = new TLegend(0.5,0.75,0.9,0.9);
  leg_hist->SetTextSize(0.033);
  leg_hist->SetHeader(legEntry);
  leg_hist -> SetFillColor(0);
  

  graph -> SetTitle("Pile-up dependence of Photon pT (data)");
  graph -> GetXaxis() -> SetTitle("#Vtx"); 
  graph -> GetYaxis() -> SetTitle("Photon pT");
  graph -> SetMaximum(230);
  graph  -> GetYaxis() -> SetTitleOffset(1.3); 
  

  TCanvas* canvas = new TCanvas("c","c",0,0,500,500);
  canvas ->cd();
  graph->Draw("AP");
  leg_hist->Draw("same");

  canvas->SaveAs("PUDependence_DATA.pdf");

  return 0;

}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int makeGraphDataBinned(){

 
  TGraphErrors * graph;
  double x[7],y[7],xErr[7],yErr[7];
  TH2D* histo2D;
  TString string;
  int pTBin = 6;

  
  string.Form("PhotonPtNVtxBinned%i_PF_data.root",pTBin);

  TFile *file = TFile::Open(string);

  if(pTBin == 0) string ="Photon Pt against Number of Vertices for 22.0 GeV < p_{T} < 36.0 GeV";
  else if (pTBin == 1) string = "Photon Pt against Number of Vertices for 36.0 GeV < p_{T} < 60.0 GeV" ;
  else if (pTBin == 2) string = "Photon Pt against Number of Vertices for 60.0 GeV < p_{T} < 88.0 GeV;" ;
  else if (pTBin == 3) string = "Photon Pt against Number of Vertices for 88.0 GeV < p_{T} < 105.0 GeV" ;
  else if (pTBin == 4) string = "Photon Pt against Number of Vertices for 105.0 GeV < p_{T} < 148.5 GeV" ;
  else if (pTBin == 5) string = "Photon Pt against Number of Vertices for 148.5 GeV < p_{T} < 165.0 GeV" ;
  else if (pTBin == 6) string = "Photon Pt against Number of Vertices for 165.0 GeV < p_{T}" ;

  file->GetObject(string,histo2D);

  string.Form("PUDependence_DATA%i.pdf",pTBin);

  TH1D * histoVtx = (TH1D*)histo2D->ProjectionY("histoVtx",0,1200);

  
  TH1D* histo1D = (TH1D*)histo2D->ProjectionX("histo1D",1,5);
  TH1D* histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histoVtxPart->SetBinContent(1,histoVtx->GetBinContent(1));
  histoVtxPart->SetBinContent(2,histoVtx->GetBinContent(2));
  histoVtxPart->SetBinContent(3,histoVtx->GetBinContent(3));
  histoVtxPart->SetBinContent(4,histoVtx->GetBinContent(4));
  histoVtxPart->SetBinContent(5,histoVtx->GetBinContent(5));
  y[0] = histo1D->GetMean();
  yErr[0] = histo1D->GetMeanError();  
  x[0] = histoVtxPart->GetMean();
  xErr[0] = histoVtxPart->GetMeanError();

  cout<<"x0 = "<<x[0]<<endl;
  cout<<"x0Err = "<<xErr[0]<<endl;

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",6,10);
  histoVtxPart->SetBinContent(6,histoVtx->GetBinContent(6));
  histoVtxPart->SetBinContent(7,histoVtx->GetBinContent(7));
  histoVtxPart->SetBinContent(8,histoVtx->GetBinContent(8));
  histoVtxPart->SetBinContent(9,histoVtx->GetBinContent(9));
  histoVtxPart->SetBinContent(10,histoVtx->GetBinContent(10));
  y[1] = histo1D->GetMean();
  yErr[1] = histo1D->GetMeanError();
  x[1] = histoVtxPart->GetMean();
  xErr[1] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",11,15);
  histoVtxPart->SetBinContent(11,histoVtx->GetBinContent(11));
  histoVtxPart->SetBinContent(12,histoVtx->GetBinContent(12));
  histoVtxPart->SetBinContent(13,histoVtx->GetBinContent(13));
  histoVtxPart->SetBinContent(14,histoVtx->GetBinContent(14));
  histoVtxPart->SetBinContent(15,histoVtx->GetBinContent(15));
  y[2] = histo1D->GetMean();
  yErr[2] = histo1D->GetMeanError();
  x[2] = histoVtxPart->GetMean();
  xErr[2] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",16,20);
  histoVtxPart->SetBinContent(16,histoVtx->GetBinContent(16));
  histoVtxPart->SetBinContent(17,histoVtx->GetBinContent(17));
  histoVtxPart->SetBinContent(18,histoVtx->GetBinContent(18));
  histoVtxPart->SetBinContent(19,histoVtx->GetBinContent(19));
  histoVtxPart->SetBinContent(20,histoVtx->GetBinContent(20));
  y[3] = histo1D->GetMean();
  yErr[3] = histo1D->GetMeanError();
  x[3] = histoVtxPart->GetMean();
  xErr[3] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",21,25);
  histoVtxPart->SetBinContent(21,histoVtx->GetBinContent(21));
  histoVtxPart->SetBinContent(22,histoVtx->GetBinContent(22));
  histoVtxPart->SetBinContent(23,histoVtx->GetBinContent(23));
  histoVtxPart->SetBinContent(24,histoVtx->GetBinContent(24));
  histoVtxPart->SetBinContent(25,histoVtx->GetBinContent(25));
  y[4] = histo1D->GetMean();
  yErr[4] = histo1D->GetMeanError();
  x[4] = histoVtxPart->GetMean();
  xErr[4] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",26,30);
  histoVtxPart->SetBinContent(26,histoVtx->GetBinContent(26));
  histoVtxPart->SetBinContent(27,histoVtx->GetBinContent(27));
  histoVtxPart->SetBinContent(28,histoVtx->GetBinContent(28));
  histoVtxPart->SetBinContent(29,histoVtx->GetBinContent(29));
  histoVtxPart->SetBinContent(30,histoVtx->GetBinContent(30));
  y[5] = histo1D->GetMean();
  yErr[5] = histo1D->GetMeanError();
  x[5] = histoVtxPart->GetMean();
  xErr[5] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",31,35);
  histoVtxPart->SetBinContent(31,histoVtx->GetBinContent(31));
  histoVtxPart->SetBinContent(32,histoVtx->GetBinContent(32));
  histoVtxPart->SetBinContent(33,histoVtx->GetBinContent(33));
  histoVtxPart->SetBinContent(34,histoVtx->GetBinContent(34));
  histoVtxPart->SetBinContent(35,histoVtx->GetBinContent(35));
  y[6] = histo1D->GetMean();
  yErr[6] = histo1D->GetMeanError();
  x[6] = histoVtxPart->GetMean();
  xErr[6] = histoVtxPart->GetMeanError();
 

  graph = new TGraphErrors(7,x,y,xErr,yErr);
  TF1* f1 = new TF1("pol1","pol1"); 
  //f1->SetParameter(0,90.);
  graph -> Fit(f1,"","",0.0,35.0);

  char legEntry[100];


  cout<<"par0 = "<<f1->GetParameter(0)<<endl;
  cout<<"par1 = "<<f1->GetParameter(1)<<endl;
  sprintf(legEntry,"%4.2f + %4.2f (#pm %4.2f) * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParError(1));

  TLegend* leg_hist       =  0 ;
  leg_hist = new TLegend(0.5,0.75,0.9,0.9);
  leg_hist->SetTextSize(0.033);
  leg_hist->SetHeader(legEntry);
  leg_hist -> SetFillColor(0);
  

  graph -> SetTitle("Pile-up dependence of Photon pT (data)");
  graph -> GetXaxis() -> SetTitle("#Vtx"); 
  graph -> GetYaxis() -> SetTitle("Photon pT");
  graph -> SetMaximum(250);
  graph  -> GetYaxis() -> SetTitleOffset(1.3); 
  

  TCanvas* canvas = new TCanvas("c","c",0,0,500,500);
  canvas ->cd();
  graph->Draw("AP");
  leg_hist->Draw("same");

  canvas->SaveAs(string);

  return 0;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int makeGraphMC(){

 
  TGraphErrors * graph;
  double x[7],y[7],xErr[7],yErr[7];
  TH2D* histo2D;
  

  TFile *file = TFile::Open("PhotonPtNVtx_2d_PF_mc.root");
  file->GetObject("Photon Pt against Number of Vertices",histo2D);

  TH1D * histoVtx = (TH1D*)histo2D->ProjectionY("histoVtx",0,1200);

  
  TH1D* histo1D = (TH1D*)histo2D->ProjectionX("histo1D",1,5);
  TH1D* histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histoVtxPart->SetBinContent(1,histoVtx->GetBinContent(1));
  histoVtxPart->SetBinContent(2,histoVtx->GetBinContent(2));
  histoVtxPart->SetBinContent(3,histoVtx->GetBinContent(3));
  histoVtxPart->SetBinContent(4,histoVtx->GetBinContent(4));
  histoVtxPart->SetBinContent(5,histoVtx->GetBinContent(5));
  y[0] = histo1D->GetMean();
  yErr[0] = histo1D->GetMeanError();  
  x[0] = histoVtxPart->GetMean();
  xErr[0] = histoVtxPart->GetMeanError();

  cout<<"x0 = "<<x[0]<<endl;
  cout<<"x0Err = "<<xErr[0]<<endl;

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",6,10);
  histoVtxPart->SetBinContent(6,histoVtx->GetBinContent(6));
  histoVtxPart->SetBinContent(7,histoVtx->GetBinContent(7));
  histoVtxPart->SetBinContent(8,histoVtx->GetBinContent(8));
  histoVtxPart->SetBinContent(9,histoVtx->GetBinContent(9));
  histoVtxPart->SetBinContent(10,histoVtx->GetBinContent(10));
  y[1] = histo1D->GetMean();
  yErr[1] = histo1D->GetMeanError();
  x[1] = histoVtxPart->GetMean();
  xErr[1] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",11,15);
  histoVtxPart->SetBinContent(11,histoVtx->GetBinContent(11));
  histoVtxPart->SetBinContent(12,histoVtx->GetBinContent(12));
  histoVtxPart->SetBinContent(13,histoVtx->GetBinContent(13));
  histoVtxPart->SetBinContent(14,histoVtx->GetBinContent(14));
  histoVtxPart->SetBinContent(15,histoVtx->GetBinContent(15));
  y[2] = histo1D->GetMean();
  yErr[2] = histo1D->GetMeanError();
  x[2] = histoVtxPart->GetMean();
  xErr[2] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",16,20);
  histoVtxPart->SetBinContent(16,histoVtx->GetBinContent(16));
  histoVtxPart->SetBinContent(17,histoVtx->GetBinContent(17));
  histoVtxPart->SetBinContent(18,histoVtx->GetBinContent(18));
  histoVtxPart->SetBinContent(19,histoVtx->GetBinContent(19));
  histoVtxPart->SetBinContent(20,histoVtx->GetBinContent(20));
  y[3] = histo1D->GetMean();
  yErr[3] = histo1D->GetMeanError();
  x[3] = histoVtxPart->GetMean();
  xErr[3] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",21,25);
  histoVtxPart->SetBinContent(21,histoVtx->GetBinContent(21));
  histoVtxPart->SetBinContent(22,histoVtx->GetBinContent(22));
  histoVtxPart->SetBinContent(23,histoVtx->GetBinContent(23));
  histoVtxPart->SetBinContent(24,histoVtx->GetBinContent(24));
  histoVtxPart->SetBinContent(25,histoVtx->GetBinContent(25));
  y[4] = histo1D->GetMean();
  yErr[4] = histo1D->GetMeanError();
  x[4] = histoVtxPart->GetMean();
  xErr[4] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",26,30);
  histoVtxPart->SetBinContent(26,histoVtx->GetBinContent(26));
  histoVtxPart->SetBinContent(27,histoVtx->GetBinContent(27));
  histoVtxPart->SetBinContent(28,histoVtx->GetBinContent(28));
  histoVtxPart->SetBinContent(29,histoVtx->GetBinContent(29));
  histoVtxPart->SetBinContent(30,histoVtx->GetBinContent(30));
  y[5] = histo1D->GetMean();
  yErr[5] = histo1D->GetMeanError();
  x[5] = histoVtxPart->GetMean();
  xErr[5] = histoVtxPart->GetMeanError();
  

  delete histoVtxPart;
  histoVtxPart = new TH1D("vtx","vtx_part",60,1,60);
  histo1D = (TH1D*)histo2D->ProjectionX("bla",31,35);
  histoVtxPart->SetBinContent(31,histoVtx->GetBinContent(31));
  histoVtxPart->SetBinContent(32,histoVtx->GetBinContent(32));
  histoVtxPart->SetBinContent(33,histoVtx->GetBinContent(33));
  histoVtxPart->SetBinContent(34,histoVtx->GetBinContent(34));
  histoVtxPart->SetBinContent(35,histoVtx->GetBinContent(35));
  y[6] = histo1D->GetMean();
  yErr[6] = histo1D->GetMeanError();
  x[6] = histoVtxPart->GetMean();
  xErr[6] = histoVtxPart->GetMeanError();
 

  graph = new TGraphErrors(7,x,y,xErr,yErr);
  TF1* f1 = new TF1("pol1","pol1"); 
  graph -> Fit(f1,"Q","",5.0,35.0);

  char legEntry[100];


  cout<<"par0 = "<<f1->GetParameter(0)<<endl;
  cout<<"par1 = "<<f1->GetParameter(1)<<endl;
  sprintf(legEntry,"%4.2f + %4.2f (#pm %4.2f) * x",f1->GetParameter(0),f1->GetParameter(1),f1->GetParError(1));

  TLegend* leg_hist       =  0 ;
  leg_hist = new TLegend(0.5,0.75,0.9,0.9);
  leg_hist->SetTextSize(0.033);
  leg_hist->SetHeader(legEntry);
  leg_hist -> SetFillColor(0);
  

  graph -> SetTitle("Pile-up dependence of Photon pT (MC)");
  graph -> GetXaxis() -> SetTitle("#Vtx"); 
  graph -> GetYaxis() -> SetTitle("Photon pT");
  graph -> SetMaximum(230);
  graph  -> GetYaxis() -> SetTitleOffset(1.3); 
  

  TCanvas* canvas = new TCanvas("c","c",0,0,500,500);
  canvas ->cd();
  graph->Draw("AP");
  leg_hist->Draw("same");

  canvas->SaveAs("PUDependence_MC.pdf");

  return 0;

}


int compareDataMC_VtxDistribution(){


  TH1D *hPhotonPtData;
  TH1D *hPhotonPtMC;
  TFile *file, *file1;  
  TH1D* hNVtxData;
  TH1D* hNVtxMC;

  char title[100];
  TH1D* hNVtxD[numTrigger];
  TH1D* hNVtxM[numTrigger];
  TCanvas* canv[numTrigger];
  TCanvas* canvas;
  for(int i=1; i<numTrigger;i++){

 
    sprintf(title,"../plots_2012/PF_L1CHS/data/root_files/VtxN%i_PFCHS_data.root",i);
    file = TFile::Open(title);  
    file->GetObject("histo",hNVtxD[i]);
    
    sprintf(title,"../plots_2012/PF_L1CHS/mc/root_files/VtxN%i_PFCHS_mc.root",i);
    file1 = TFile::Open(title);
    file1->GetObject("histo",hNVtxM[i]);
    
    hNVtxM[i]->Scale(hNVtxD[i]->Integral()/hNVtxM[i]->Integral()); 
    
    sprintf(title,"Number of Vertices (data and MC) pt %i  NEW PU Dist",i);
    canv[i] = DrawComparison(hNVtxD[i], hNVtxM[i],title,"","#Vtx", 1);
    canv[i] ->SetLogy();
    sprintf(title,"plots/NVtxComparison%i.pdf",i);
    canv[i] ->Print(title);
  
  }

  delete file;
  delete file1;
  
  return 0;
}

int makeRhoGraph(){

  char title[100];
  TH1D* hRhoD[numTrigger];
  TH1D* hRhoM[numTrigger];
  TCanvas* canv[numTrigger];
  TFile *file;
  TFile *file1;
  for(int i=0; i<numTrigger;i++){

    sprintf(title,"Rho%i_PF_data.root",i);
    file = TFile::Open(title);
  
    sprintf(title,"hRhoPtBinned%i",i);
    file->GetObject(title,hRhoD[i]);
    
    sprintf(title,"Rho%i_PF_mc.root",i);
    file1 = TFile::Open(title);
    
    sprintf(title,"hRhoPtBinned%i",i);
    file1->GetObject(title,hRhoM[i]);
    
 
    hRhoM[i]->Scale(hRhoD[i]->Integral()/hRhoM[i]->Integral()); 
    
    sprintf(title,"Rho (data and MC) pt %i ",i);
    canv[i] = DrawComparison(hRhoD[i], hRhoM[i],title,"","Rho", 1);
    
    sprintf(title,"RhoComparison%i.pdf",i);
    
    canv[i]->Print(title);
   
    delete file;
    delete file1;
    
  }

  return 0;
}

int makeOtherRhoGraph(){

  char title[100];
  TH1D* hRhoD[40];
  TH1D* hRhoM[40];
  TCanvas* canv[40];
  TFile *file;
  TFile *file1;
  for(int i=1; i<40;i++){
    
    sprintf(title,"RhoVtxBinned/RhoVtxBinned%i_PF_data.root",i);
    file = TFile::Open(title);
   
    sprintf(title,"hRhoVtxBinned%i",i);
    file->GetObject(title,hRhoD[i]);
     
    sprintf(title,"RhoVtxBinned/RhoVtxBinned%i_PF_mc.root",i);
    file1 = TFile::Open(title);
     
    sprintf(title,"hRhoVtxBinned%i",i);
    file1->GetObject(title,hRhoM[i]);
  
 
    hRhoM[i]->Scale(hRhoD[i]->Integral()/hRhoM[i]->Integral()); 
     
    sprintf(title,"Rho (data and MC) for #Vtx = %i ",i);
    canv[i] = DrawComparison(hRhoD[i], hRhoM[i],title,"","Rho", 1);
     
    sprintf(title,"RhoVtxBinned/RhoComparisonVtxBinned%i.pdf",i);
      
    canv[i]->Print(title);
     
    delete file;
    delete file1;
    
  }

  return 0;
}

int makeTriggerEffGraph(){

  char title[100];
  TH1D* hTriggerBefore[numTrigger];
  TH1D* hTriggerAfter[numTrigger];
  TH1D* hTrigger[numTrigger];
  TCanvas* canv[numTrigger];
  TFile *file;
  TFile *file1;
  for(int i=0; i<numTrigger;i++){


    sprintf(title,"TriggerEffBefore%i_PF_mc.root",i);
    file = TFile::Open(title);
    
    sprintf(title,"hTriggerEffPtBinnedBefore%i",i);
    file->GetObject(title,hTriggerBefore[i]);
     
    sprintf(title,"TriggerEffAfter%i_PF_mc.root",i);
    file1 = TFile::Open(title);
    
    sprintf(title,"hTriggerEffPtBinnedAfter%i",i);
    file1->GetObject(title,hTriggerAfter[i]);
    
    hTrigger[i] = new TH1D("ratio","ratio",60,0,60);

    hTrigger[i]->Divide(hTriggerAfter[i],hTriggerBefore[i]);
    canv[i] = new TCanvas("ratio","ratio",500,500);
    canv[i] -> cd();

    hTrigger[i]->SetTitle("Trigger efficiency in MC");

    hTrigger[i]->SetXTitle("#Vtx");

    hTrigger[i]->SetYTitle("#events pass Trigger / # all events");

    hTrigger[i]->Draw();
   
    sprintf(title,"TriggerEff%i.pdf",i);
    canv[i]->Print(title);
    
    delete file;
    delete file1;
    
  }

  return 0;
}


TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool isData)
{
   double MinX = selection->GetXaxis()->GetXmin();
   double MaxX = selection->GetXaxis()->GetXmax();
   double MaxY = selection->GetMaximum();
   double YRangeMax = 2*MaxY;
   TString titlePrediction;
   TString titleSelection;
   TString RatioTitle;
   
   if( isData ){
      titlePrediction = "Data";
      titleSelection = "MC";
      RatioTitle = "(Data-MC)/MC";
   }
   else {
      titlePrediction = "Data-driven Pred. from MC";
      titleSelection = "MC Expectation";
      RatioTitle = "(Pred-MC)/MC";
   }

   //static Int_t c_LightBrown   = TColor::GetColor( "#D9D9CC" );
   static Int_t c_LightGray    = TColor::GetColor( "#DDDDDD" );

   prediction->SetAxisRange(MinX, MaxX, "X");
   prediction->GetYaxis()->SetRangeUser(0.005, YRangeMax);
   prediction->SetMarkerStyle(20);
   prediction->SetMarkerSize(0.9);
   prediction->SetMarkerColor(kBlack);
   prediction->SetXTitle(xTitle);
   prediction->SetYTitle("Events");
   selection->SetAxisRange(MinX, MaxX, "X");
   selection->GetYaxis()->SetRangeUser(0.05, YRangeMax);
   // selection->SetFillColor(c_LightBrown);
   selection->SetFillColor(c_LightGray);
   selection->SetTitle("");
   selection->SetXTitle(xTitle);
   selection->SetYTitle("Events");
   TCanvas *c = new TCanvas("ca", "Comparison and ratio of two histos", 700, 700);
   TPad *pad1 = new TPad("pad1a", "pad1a", 0, 0.35, 1, 1);
   pad1->SetLogy();
   pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();
  
   selection->DrawCopy("hist");
   prediction->Draw("same");
   selection->SetFillColor(kAzure-3);
   selection->SetFillStyle(3354);
   selection->DrawCopy("e2same");
   selection->SetFillStyle(1001);
   //  selection->SetFillColor(c_LightBrown);
   selection->SetFillColor(c_LightGray);

   TLegend* leg1 = new TLegend(0.48, 0.63, 0.95, 0.83);
   leg1->SetFillStyle(0);
   leg1->SetLineStyle(1);
   leg1->SetTextFont(42);
   leg1->SetTextSize(0.04);
   leg1->AddEntry(selection, titleSelection, "lf");
   leg1->AddEntry(prediction, titlePrediction, "lep");
   leg1->Draw("same");
 
   TPaveText* pt = new TPaveText(0.11, 0.98, 0.95, 0.86, "NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.045);
   pt->AddText(Title);
   pt->AddText(LumiTitle);
   pt->Draw();
   c->cd();
   TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 0.35);
   pad2->SetTopMargin(0);
   pad2->Draw();
   pad2->cd();
   TH1D* r = new TH1D(*selection);
   r->SetTitle("");
   r->SetLabelSize(0.08, "XYZ");
   r->SetLabelOffset(0.01, "XYZ");
   r->SetTitleSize(0.09, "XYZ");
   r->SetTitleOffset(0.65, "Y");
   r->SetTickLength(0.05);
   r->SetYTitle(RatioTitle);
   r->SetStats(0);
   r->SetMarkerStyle(20);
   r->SetMarkerSize(0.9);
   r->SetMarkerColor(kBlack);
   r->Reset();
   r->Add(prediction, 1);
   r->Add(selection, -1);
   r->Divide(selection);
   r->SetMaximum(1.2);
   r->SetMinimum(-1.2);
   //r->SetMaximum(0.5);
   //r->SetMinimum(-0.5);
   r->Draw("ep");
   TLine l;
   l.DrawLine(MinX, 0., MaxX, 0.);
   c->cd();
   return c;
}
