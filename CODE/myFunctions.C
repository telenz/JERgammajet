#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TStyle.h"
#include "myFunctions.h"
#include "myDeclarations.h"

int histoInt = 0;

void plotTH1(TH1* histo,const TString &filename, bool logy){

  TString tot_filename;
  TFile* f;

  TCanvas* canvas = new TCanvas(filename,filename,0,0,500,500);
  
  canvas -> SetLeftMargin(0.17);
 
  histo ->GetYaxis()->SetTitleOffset(1.5); 
  //histo -> GetYaxis()->SetNoExponent(1);
  canvas -> cd();
  histo -> Draw(); 
  if(logy) canvas ->SetLogy();
  
  tot_filename = PDFPath + filename + DataType + ".pdf";
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  //canvas -> SaveAs(tot_filename);
  gErrorIgnoreLevel = oldLevel;
  
  tot_filename = RootPath + filename + DataType + ".root";
  f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(histo,"histo");
  f ->Close();
  delete f;
  delete canvas;
  
}

void plotTH2(TH2* histo,const TString &filename){
  
  TString tot_filename;
  TFile* f; 
  TCanvas* canvas = new TCanvas(filename,filename,0,0,500,500);
  canvas -> SetLeftMargin(0.17);
  canvas -> SetRightMargin(0.17);
  canvas ->SetLogz();
  histo -> GetYaxis()->SetTitleOffset(1.5); 
  histo -> GetZaxis()->CenterTitle();
  //histo -> GetYaxis()->SetDecimals();
  canvas -> cd();
  histo -> Draw("COLZ");
  tot_filename = PDFPath + filename + DataType + ".pdf";
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  canvas -> SaveAs(tot_filename);
  gErrorIgnoreLevel = oldLevel;

  tot_filename = RootPath + filename + DataType + ".root";
  f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(histo,"histo");
  f ->Close();
  delete f;
  delete canvas;
  }

void plotTGraphErrors(TGraphErrors* graph, const TString &filename, const TString &graphTitle,const TString &xTitle,const TString &yTitle, double x_low, double x_up, double y_low, double y_up,TString legEntry, bool legend){
   
  TString tot_filename;
  TFile* f;

  TCanvas* canvas = new TCanvas(filename,filename,0,0,400,400);
  canvas -> SetLeftMargin(0.17);
  graph  -> SetTitle(graphTitle);
  graph  -> GetXaxis() -> SetTitle(xTitle); 
  graph  -> GetYaxis() -> SetTitle(yTitle);
  graph  -> GetYaxis() -> SetTitleOffset(1.5); 
  graph  -> GetXaxis() -> SetTitleOffset(1.1); 
  graph  -> GetYaxis() -> SetRangeUser(y_low,y_up); 
  graph  -> GetXaxis() -> SetLimits(x_low,x_up); 


  canvas -> cd();
  graph  -> Draw("AP"); 

  if(legend){
  TLegend* leg_hist       =  0 ;
  leg_hist = new TLegend(0.4,0.75,0.9,0.9);
  leg_hist->SetTextSize(0.033);
  leg_hist->SetHeader(legEntry);
  leg_hist -> SetFillColor(0);
  leg_hist->Draw("same");
  }

  tot_filename = PDFPath + filename + DataType + MethodType +".pdf";
  //tot_filename = "test.pdf";
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  canvas -> SaveAs(tot_filename);
  gErrorIgnoreLevel = oldLevel;
  tot_filename = RootPath + filename + DataType + MethodType +".root";
  //tot_filename = "test.root";
  f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(graph,"Graph");
  f -> Close();
  delete f;
  delete canvas;

}

void saveObject(TObject* object, TString filename, TString objectName){

  TFile* f;
  f = new TFile(filename,"RECREATE");
  f -> WriteTObject(object,objectName);
  f->Close();
  delete f;
  
}


TH1D* createTH1(TH1D *histo,const TString &HistoName,const TString &HistoTitle,int nBin, double bound_low, double bound_up,const TString &xTitle){


  TString HistoName2;
  if(HistoName.EqualTo("")){
    HistoName2.Form("histo%i",histoInt);
    histoInt = histoInt + 1;
  }
  else HistoName2 = HistoName;
  

  histo = new TH1D(HistoName2,HistoTitle,nBin,bound_low,bound_up);
  histo -> Sumw2();
  histo -> SetXTitle(xTitle);
    
  

  return histo;
}

TH2D* createTH2(TH2D *histo,const TString &HistoTitle,int XnBin,double Xbound_low,double Xbound_up,int YnBin,double Ybound_low,double Ybound_up,const TString &xTitle,const TString &yTitle){
 
      TString histoName;
      std::stringstream out; 
 
      out<< histoInt;
      histoName = out.str();

      histo = new TH2D(HistoTitle,HistoTitle, XnBin, Xbound_low, Xbound_up, YnBin, Ybound_low, Ybound_up);
      histo -> Sumw2();
      histo -> SetXTitle(xTitle);
      histo -> SetYTitle(yTitle);
      //histo -> SetZTitle("# events");
      
      histoInt = histoInt + 1;

  return histo;
}


TH1D* readTH1(const TString &fileName, const TString &histName, const TString &newHistName) {
    TFile f(fileName,"READ");
    TH1D *h = 0;
    f.GetObject(histName,h);
    if( h ) {
      h->SetDirectory(0);
   
      if( newHistName.Length() ) h->SetName(newHistName);
    } else {
      std::cerr << "ERROR in FileOps::readTH1: Histogram with name '" << histName << "' does not exist in file '" << fileName << "'\n.";
      f.Close();
      exit(-1);
    }
    f.Close();
    
    return h;
}




TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName) {
    TFile f(fileName,"READ");
    TGraphErrors *g = 0;
    f.GetObject(gName,g);
    if( g ) {
      if( newGName.Length() ) g->SetName(newGName);
    } else {
      std::cerr << "ERROR in FileOps::readTGraph: TGraphErrors with name '" << gName << "' does not exist in file '" << fileName << "'\n.";
      f.Close();
      exit(-1);
    }
    f.Close();
    
    return g;
    }
