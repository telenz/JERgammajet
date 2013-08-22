
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
#include "TColorWheel.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TROOT.h"
#include "../CODE/myDeclarations.h"
#include "plotStyle.h"

TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool isData);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotFlavorComposition(){
  

  //TString definition = "algo";
  TString definition = "phys";
  
  TString pwdName = "../plots_2012/PF_L1CHS/mc/root_files/";

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetMarkerStyle(20);

  TH2D *charm2D, *bottom2D, *gluon2D, *lightQuark2D, *nonDefined2D;
  TH1D *charm[2], *bottom[2], *gluon[2], *lightQuark[2], *together[2], *nonDefined[2], *allQuarks[2], *togetherWoUndefined[2];
  TH1D *charmComp[2], *bottomComp[2], *gluonComp[2], *lightQuarkComp[2], *nonDefinedComp[2], *allQuarksComp[2], *gluonWoUndefinedComp[2];
  TString fileName;  
  TFile *file;
  TCanvas* canvas[2] = {0};  
  TLatex*  info[2]   = {0};
  TLegend* legend[2] = {0};
  gStyle->SetOptStat("");
  
   
  // Read Files
  fileName =  pwdName + "hPhotonPtJetEtaFlavorBottom_" + definition + "_PFCHS_mc.root";
  file     =  new TFile(fileName);
  bottom2D =  (TH2D*) gDirectory->Get("histo");
  bottom2D -> SetDirectory(0);
  bottom[0]  = (TH1D*) bottom2D->ProjectionX("bottom0",0,bottom2D->GetYaxis()->FindBin(1.3));
  bottom[1]  = (TH1D*) bottom2D->ProjectionX("bottom1",bottom2D->GetYaxis()->FindBin(1.3),-1);
  cout<<"bottom2D->GetYaxis()->FindBin(1.3) = "<<bottom2D->GetYaxis()->FindBin(1.3)<<endl;
  cout<<"bottom2D->GetYaxis()->FindBin(5.0) = "<<bottom2D->GetYaxis()->FindBin(5.0)<<endl;
  for(int i= 0; i<2; i++) bottom[i] -> Rebin(20);
  for(int i= 0; i<2; i++) bottom[i] -> SetDirectory(0);
  delete file;

  fileName =  pwdName + "hPhotonPtJetEtaFlavorCharm_" + definition + "_PFCHS_mc.root";
  file     =  new TFile(fileName);
  charm2D  =  (TH2D*) gDirectory->Get("histo");
  charm2D  -> SetDirectory(0);
  charm[0] = charm2D->ProjectionX("charm0",0,charm2D->GetYaxis()->FindBin(1.3),"");  
  charm[1] = charm2D->ProjectionX("charm1",charm2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++) charm[i]  ->Rebin(20);
  for(int i= 0; i<2; i++) charm[i]  -> SetDirectory(0);
  delete file;

  fileName =  pwdName + "hPhotonPtJetEtaFlavorGluon_" + definition + "_PFCHS_mc.root";
  file     =  new TFile(fileName);
  gluon2D  =  (TH2D*) gDirectory->Get("histo");
  gluon2D  -> SetDirectory(0);
  gluon[0]    = gluon2D->ProjectionX("gluon0",0,gluon2D->GetYaxis()->FindBin(1.3),"");
  gluon[1]    = gluon2D->ProjectionX("gluon1",gluon2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++)  gluon[i]    -> Rebin(20);
  for(int i= 0; i<2; i++)  gluon[i]      -> SetDirectory(0); 
  delete file;
  
  fileName     =  pwdName + "hPhotonPtJetEtaFlavorLightQuarks_" + definition + "_PFCHS_mc.root";
  file         =  new TFile(fileName);
  lightQuark2D =  (TH2D*) gDirectory->Get("histo");
  lightQuark2D -> SetDirectory(0);
  lightQuark[0] = lightQuark2D->ProjectionX("lightQuark0",0,lightQuark2D->GetYaxis()->FindBin(1.3),"");
  lightQuark[1] = lightQuark2D->ProjectionX("lightQuark1",lightQuark2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++)   lightQuark[i] -> Rebin(20);
  for(int i= 0; i<2; i++)   lightQuark[i] -> SetDirectory(0);


  fileName     =  pwdName + "hPhotonPtJetEtaFlavorNonDefined_" + definition + "_PFCHS_mc.root";
  file         =  new TFile(fileName);
  nonDefined2D =  (TH2D*) gDirectory->Get("histo");
  nonDefined2D -> SetDirectory(0);
  nonDefined[0] = nonDefined2D->ProjectionX("nonDefined0",0,nonDefined2D->GetYaxis()->FindBin(1.3),"");
  nonDefined[1] = nonDefined2D->ProjectionX("nonDefined1",nonDefined2D->GetYaxis()->FindBin(1.3),-1,"");
  for(int i= 0; i<2; i++)  nonDefined[i] -> Rebin(20);
  for(int i= 0; i<2; i++)  nonDefined[i] -> SetDirectory(0);
  delete file;
     
  for(int i= 0; i<2; i++){

    allQuarks[i] = (TH1D*) bottom[i]->Clone(Form("allQuarks%i",i));
    allQuarks[i] -> Add(charm[i]);
    allQuarks[i] -> Add(lightQuark[i]);

    together[i] = (TH1D*) bottom[i]->Clone(Form("together%i",i));
    together[i] -> Add(charm[i]);
    together[i] -> Add(gluon[i]);
    together[i] -> Add(lightQuark[i]);
    together[i] -> Add(nonDefined[i]);

    togetherWoUndefined[i] = (TH1D*) bottom[i]->Clone(Form("together%i",i));
    togetherWoUndefined[i] -> Add(charm[i]);
    togetherWoUndefined[i] -> Add(gluon[i]);
    togetherWoUndefined[i] -> Add(lightQuark[i]);

    togetherWoUndefined[i]->Rebin(togetherWoUndefined[i]->GetNbinsX());
    allQuarks[i]->Rebin(allQuarks[i]->GetNbinsX());

    cout<<"number of bins together = "<<togetherWoUndefined[i]->GetNbinsX()<<endl;
    cout<<"number of bins allQuarks = "<<allQuarks[i]->GetNbinsX()<<endl;
    
    bottomComp[i]     = (TH1D*) bottom[i]->Clone(Form("bottom%i",i+3));
    bottomComp[i]     -> Divide(together[i]);
    charmComp[i]      = (TH1D*) charm[i]->Clone(Form("charm%i",i+3));
    charmComp[i]      -> Divide(together[i]);
    gluonComp[i]      = (TH1D*) gluon[i]->Clone(Form("gluon%i",i+3));
    gluonComp[i]      -> Divide(together[i]);
    lightQuarkComp[i] = (TH1D*) lightQuark[i]->Clone(Form("lightQuark%i",i+3));
    lightQuarkComp[i] -> Divide(together[i]);
    nonDefinedComp[i] = (TH1D*) nonDefined[i]->Clone(Form("nonDefined%i",i+3));  
    nonDefinedComp[i] -> Divide(together[i]);
    allQuarksComp[i]  = (TH1D*) allQuarks[i]->Clone(Form("allQuarks%i",i+3));
    allQuarksComp[i] -> Divide(togetherWoUndefined[i]);
    gluon[i]->Rebin(gluon[i]->GetNbinsX());
    gluonWoUndefinedComp[i] = (TH1D*) gluon[i]->Clone(Form("gluon%i",i+3));
    gluonWoUndefinedComp[i] -> Divide(togetherWoUndefined[i]);

    cout<<endl<<endl<<"Gluon Flavor Fraction for systematics  = "<<gluonWoUndefinedComp[i]->GetBinContent(1)<<endl;
    cout<<"Quark Flavor Fraction for systematics  = "<<allQuarksComp[i]->GetBinContent(1)<<endl;
   
    // Set Marker Style 
    lightQuarkComp[i] -> SetMarkerStyle(20);
    gluonComp[i]      -> SetMarkerStyle(24);
    charmComp[i]      -> SetMarkerStyle(21);
    bottomComp[i]     -> SetMarkerStyle(23);
    nonDefinedComp[i] -> SetMarkerStyle(22);
    lightQuarkComp[i] -> SetMarkerColor(2);
    gluonComp[i]      -> SetMarkerColor(1);
    charmComp[i]      -> SetMarkerColor(3);
    bottomComp[i]     -> SetMarkerColor(4);    
    nonDefinedComp[i] -> SetMarkerColor(5);
    lightQuarkComp[i] -> SetLineColor(2);
    gluonComp[i]      -> SetLineColor(1);
    charmComp[i]      -> SetLineColor(3);
    bottomComp[i]     -> SetLineColor(4);
    nonDefinedComp[i] -> SetLineColor(5);

    canvas[i] = new TCanvas(Form("canvas%i",i),"canvas",0,0,500,500);
    canvas[i] ->cd();
    //canvas[i] ->SetLogy();

    // Plot the results
    lightQuarkComp[i] -> SetTitle("Flavor Composition in #gamma + Jet Events");
    lightQuarkComp[i] -> GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    lightQuarkComp[i] -> SetTitleSize(0.06,"X");
    lightQuarkComp[i] -> GetYaxis()->SetTitle("Flavor Fraction");
    lightQuarkComp[i] -> SetTitleSize(0.06,"Y");
    lightQuarkComp[i] -> GetYaxis() -> SetTitleOffset(1.1);
    lightQuarkComp[i] -> SetMinimum(0.0);
    lightQuarkComp[i] -> SetMaximum(1.2);
    lightQuarkComp[i] -> Draw(); 
    //quarksComp -> Draw(); 
    gluonComp[i]      -> Draw("same");
    charmComp[i]      -> Draw("same");
    bottomComp[i]     -> Draw("same");
    nonDefinedComp[i] -> Draw("same");

    legend[i] = new TLegend(0.5,0.65,0.9,0.9);
    legend[i] -> SetTextSize(0.043);
    legend[i] -> SetFillColor(0);
    
    legend[i] -> AddEntry(lightQuarkComp[i],"uds quarks", "pfl");
    legend[i] -> AddEntry(charmComp[i],"charm quarks", "pfl");
    legend[i] -> AddEntry(bottomComp[i],"bottom quarks", "pfl");
    legend[i] -> AddEntry(gluonComp[i],"gluons", "pfl");
    legend[i] -> AddEntry(nonDefinedComp[i],"not defined", "pfl");
    legend[i] -> Draw("same");


    info[i]   = new TLatex();
    info[i] -> SetTextFont(132);
    info[i] -> SetNDC();
    info[i] -> SetTextSize(0.043);
  
    if(i==0)      info[i] -> DrawLatex(0.2,0.75,  "#bf{0.0 < |#eta^{Jet}| < 1.3}");
    else if(i==1) info[i] -> DrawLatex(0.2,0.75,  "#bf{1.3 < |#eta^{Jet}| < 2.3}");

    canvas[i] -> SetBottomMargin(0.15);
    canvas[i] -> SetLeftMargin(0.15);  
    if(i==0) canvas[i] -> SaveAs(Form("flavorFraction_barrel_" + definition + ".pdf",i+1));
    else if(i==1) canvas[i] -> SaveAs(Form("flavorFraction_endcap_" + definition + ".pdf",i+1));


  }
  //quarksComp = (TH1D*) quarks->Clone("histo");
  //quarksComp -> Add(charm);
  //quarksComp -> Add(bottom);
  //quarksComp -> Divide(together); 

  return 0;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotResponse(){
  
  gStyle -> SetTitleSize(0.05,"X"); 
  gStyle -> SetTitleSize(0.05,"Y"); 
  gStyle -> SetLabelSize(0.04,"X");
  gStyle -> SetLabelSize(0.04,"Y");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.2,"X");
  gStyle -> SetTitleOffset(1.2,"Y");
  gStyle -> SetOptLogy();
  gStyle -> SetHistLineWidth(1);
  gROOT->ForceStyle(); 
  gStyle->SetOptStat("");

  TH1D *histoPhoton1, *histoJet1, *histoPhoton2 , *histoJet2, *histoJet3, *histoTogether1, *histoTogether2;   
  TCanvas* canvas1, *canvas2, *canvas3, *canvas4, *canvas5;
  TFile *file;
  TString fileName;
  TLatex*  info;
      
  
  int etaB = 1;
  int alphaB = 6;
  int ptB = 12;
  
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Full Resolution

  canvas1 = new TCanvas("canvas1","canvas",0,0,500,500);
  canvas1 ->cd();

  fileName = Form("root_files/response_photon_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoPhoton1);
  histoPhoton1 -> SetDirectory(0); 
  delete file; 

  fileName = Form("root_files/response_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoJet1);
  histoJet1    -> SetDirectory(0); 
  delete file; 

  fileName = Form("root_files/response_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoTogether1);
  histoTogether1    -> SetDirectory(0); 
  delete file; 


  histoTogether1 -> Add(histoPhoton1);
  histoTogether1 -> Scale(histoTogether1->GetEntries()/histoTogether1->Integral());
  histoTogether1 -> Rebin(50);

  histoTogether1 -> GetXaxis()->SetTitle("p_{T}^{reco. jet}/p_{T}^{#gamma}");
  histoTogether1 -> GetYaxis()->SetTitle("# Events");
  histoTogether1 -> SetTitle("Measured Response");
  histoTogether1 -> SetMinimum(1.);
  histoTogether1 -> SetMaximum(histoTogether1->GetMaximum()*20.);

  cout<<"histoTogether->GetMarkerColor() = "<<histoTogether1->GetMarkerColor()<<endl;
  
  histoTogether1 -> SetMarkerColor(1);
  histoTogether1 -> SetLineColor(1);

  histoTogether1 -> Draw(); 


  // Draw info boxes
  info   = new TLatex();
  info -> SetTextFont(132);
  info -> SetNDC();
  info -> SetTextSize(0.041);
  
  if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
  
  canvas1->SaveAs("plots/fullResponseExample.pdf");
  //-------------------------------------------------------------------------------------------------------------------------

  canvas4 = new TCanvas("canvas4","canvas",0,600,500,500);
  canvas4 ->cd();

  histoTogether1 -> Draw(); 
  
  double scaleFactor = histoJet1->Integral() + histoPhoton1->Integral(); 

  histoJet1 -> Scale(histoTogether1->GetEntries()/scaleFactor);

  histoJet1 -> Rebin(50);
  histoJet1 -> SetMarkerColor(8);
  histoJet1 -> SetLineColor(8);
  histoJet1 -> Draw("same");

  histoPhoton1 -> Scale(histoTogether1->GetEntries()/scaleFactor);
  histoPhoton1 -> Rebin(50);
  histoPhoton1 -> SetMarkerColor(2);
  histoPhoton1 -> SetLineColor(2);

  histoPhoton1 -> Draw("same");


  // Draw info boxes
  info   = new TLatex();
  info -> SetTextFont(132);
  info -> SetNDC();
  info -> SetTextSize(0.041);
  
  //if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  //else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  //info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  //info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));


  TLegend *legend1 = new TLegend(0.25,0.75,0.9,0.9);
  legend1 -> SetTextSize(0.033);
  legend1 -> SetFillColor(0);
  legend1 -> SetLineWidth(2);
  legend1 -> AddEntry(histoJet1, "2nd jet in leading jet hemisphere");
  legend1 -> AddEntry(histoPhoton1, "2nd jet in photon hemisphere");
  legend1 -> Draw("same");

  canvas4->SaveAs("plots/fullResponseAndContributionsExample.pdf");

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Imbalance
  canvas2 = new TCanvas("canvas2","canvas",500,0,500,500);
  canvas2 ->cd();
  
  fileName = Form("root_files/response_imbalance_photon_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoPhoton2);
  histoPhoton2 -> SetDirectory(0); 
  delete file; 

  fileName = Form("root_files/response_imbalance_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoJet2);
  histoJet2    -> SetDirectory(0); 
  delete file; 


  fileName = Form("root_files/response_imbalance_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoTogether2);
  histoTogether2    -> SetDirectory(0); 
  delete file; 


  histoTogether2 -> Add(histoPhoton2);
  histoTogether2 -> Scale(histoTogether2->GetEntries()/histoTogether2->Integral());
  histoTogether2 -> Rebin(50);

  histoTogether2 -> GetXaxis()->SetTitle("p_{T}^{gen. jet}/p_{T}^{#gamma}");
  histoTogether2 -> GetYaxis()->SetTitle("# Events");
  histoTogether2 -> SetTitle("Imbalance");
  histoTogether2 -> SetMinimum(1.);
  histoTogether2 -> SetMaximum(histoTogether2->GetMaximum()*20.);

  histoTogether2 -> SetMarkerColor(4);
  histoTogether2 -> SetLineColor(4);
  
  histoTogether2 -> Draw(); 


  // Draw info boxes
  info   = new TLatex();
  info -> SetTextFont(132);
  info -> SetNDC();
  info -> SetTextSize(0.041);
  
  if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
  
  canvas2->SaveAs("plots/imbalanceExample.pdf");

  //-------------------------------------------------------------------------------------------------------------------------

  canvas5 = new TCanvas("canvas5","canvas",500,600,500,500);
  canvas5 ->cd();

  histoTogether2 -> Draw(); 
  
  scaleFactor = histoJet2->Integral() + histoPhoton2->Integral(); 

  histoJet2 -> Scale(histoTogether2->GetEntries()/scaleFactor);

  histoJet2 -> Rebin(50);
  histoJet2 -> SetMarkerColor(8);
  histoJet2 -> SetLineColor(8);
  histoJet2 -> Draw("same");

  histoPhoton2 -> Scale(histoTogether2->GetEntries()/scaleFactor);
  histoPhoton2 -> Rebin(50);
  histoPhoton2 -> SetMarkerColor(2);
  histoPhoton2 -> SetLineColor(2);
  histoPhoton2 -> Draw("same");


  // Draw info boxes
  info   = new TLatex();
  info -> SetTextFont(132);
  info -> SetNDC();
  info -> SetTextSize(0.041);
  
  //if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  //else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  //info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  //info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));

  TLegend *legend2 = new TLegend(0.25,0.75,0.9,0.9);
  legend2 -> SetTextSize(0.033);
  legend2 -> SetFillColor(0);
  legend2 -> SetLineWidth(2);
  legend2 -> AddEntry(histoJet2, "2nd jet in leading jet hemisphere");
  legend2 -> AddEntry(histoPhoton2, "2nd jet in photon hemisphere");
  legend2 -> Draw("same");
  
  canvas5->SaveAs("plots/imbalanceAndContributionsExample.pdf");

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Inrinsic
  canvas3 = new TCanvas("canvas3","canvas",1000,0,500,500);
  canvas3 ->cd();
  
  fileName = Form("root_files/response_intrinsic_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  file = TFile::Open(fileName);
  file        -> GetObject("histo",histoJet3);
  histoJet3 -> SetDirectory(0); 
  delete file; 

  histoJet3 ->Scale(histoJet3->GetEntries()/histoJet3->Integral());
  histoJet3 -> Rebin(50);

  histoJet3 -> GetXaxis()->SetTitle("p_{T}^{recon. jet}/p_{T}^{gen. jet}");
  histoJet3 -> GetYaxis()->SetTitle("# Events");
  histoJet3 -> SetTitle("Intrinsic");
  histoJet3 -> SetMinimum(1.);
  histoJet3 -> SetMaximum(histoJet3->GetMaximum()*20.);

  histoJet3 -> SetMarkerColor(46);
  histoJet3 -> SetLineColor(46);
  
  histoJet3 -> Draw(); 


  // Draw info boxes
  info   = new TLatex();
  info -> SetTextFont(132);
  info -> SetNDC();
  info -> SetTextSize(0.041);
  
  if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
  
  canvas3->SaveAs("plots/intrinsicExample.pdf");
  
  return 0;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotAllResponseHistograms_MC(){
  
  gStyle -> SetTitleSize(0.05,"X"); 
  gStyle -> SetTitleSize(0.05,"Y"); 
  gStyle -> SetLabelSize(0.04,"X");
  gStyle -> SetLabelSize(0.04,"Y");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.2,"X");
  gStyle -> SetTitleOffset(1.2,"Y");
  gStyle -> SetOptLogy();
  gStyle -> SetHistLineWidth(1);
  gROOT->ForceStyle(); 
  gStyle->SetOptStat("");

  const int etaB = 1;
  const int alphaB = 6;
  const int ptB = 12;

  TH1D *histoJet1;
  //TH1D *histoPhoton1;   
  TCanvas *canvas1[ptB][etaB][alphaB][2];
  TFile *file;
  TString fileName;
  TLatex*  info;
      
  
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Full Resolution

  for(int j=0; j<etaB; j++){
    for(int i=0; i<ptB; i++){
      for(int k=0; k<alphaB; k++){
	/*
	file = 0;
	fileName = Form("root_files/response_intrinsic_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",i+1,j+1,k+1);
	//cout<<"filename = "<<fileName<<endl;
	file = TFile::Open(fileName);
	if(file == 0) continue;
	file        -> GetObject("histo",histoPhoton1);
	histoPhoton1 -> SetDirectory(0);
	delete file; 
	if(histoPhoton1->GetEntries() < 100) continue;
	*/
	
	file = 0;
	fileName = Form("root_files/response_intrinsic_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",i+1,j+1,k+1);
	//cout<<"filename = "<<fileName<<endl;
	file = TFile::Open(fileName);
	if(file == 0) continue;
	file        -> GetObject("histo",histoJet1);
	histoJet1    -> SetDirectory(0); 
	delete file; 
	if(histoJet1->GetEntries() < 100) continue;
	

	canvas1[i][j][k][0] = new TCanvas(Form("canvas0%i%i%i",i,j,k),"canvas",0,0,500,500);
	canvas1[i][j][k][0] -> cd();
	histoJet1 -> Scale(histoJet1->GetEntries()/histoJet1->Integral());
	histoJet1 -> Rebin(50);
	histoJet1 -> GetXaxis()->SetTitle("p_{T}^{recon. jet}/p_{T}^{#gamma}");
	histoJet1 -> GetYaxis()->SetTitle("# Events");
	histoJet1 -> SetTitle("Measured Response - Jet Hemisphere");
	histoJet1 -> SetMinimum(1.);
	histoJet1 -> SetMaximum(histoJet1->GetMaximum()*20.);	
	histoJet1 -> Draw(); 
	// Draw info boxes
	info   = new TLatex();
	info -> SetTextFont(132);
	info -> SetNDC();
	info -> SetTextSize(0.041);
	
	if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
	else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
	info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
	info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
	
	fileName = Form("plotsResponse/response_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.pdf",i,j,k);
	canvas1[i][j][k][0]->SaveAs(fileName);

	/*
	canvas1[i][j][k][1] = new TCanvas(Form("canvas1%i%i%i",i,j,k),"canvas",0,0,500,500);
	canvas1[i][j][k][1] -> cd();
	histoPhoton1 -> Scale(histoPhoton1->GetEntries()/histoPhoton1->Integral());
	histoPhoton1 -> Rebin(50);
	histoPhoton1 -> GetXaxis()->SetTitle("p_{T}^{recon. jet}/p_{T}^{#gamma}");
	histoPhoton1 -> GetYaxis()->SetTitle("# Events");
	histoPhoton1 -> SetTitle("Measured Response - Photon Hemisphere");
	histoPhoton1 -> SetMinimum(1.);
	histoPhoton1 -> SetMaximum(histoPhoton1->GetMaximum()*20.);	
	histoPhoton1 -> Draw(); 

	// Draw info boxes
	info   = new TLatex();
	info -> SetTextFont(132);
	info -> SetNDC();
	info -> SetTextSize(0.041);
	
	if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
	else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
	info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
	info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
	
	fileName = Form("plotsResponse/response_photon_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.pdf",i+1,j+1,k+1);
	canvas1[i][j][k][1]->SaveAs(fileName);
	*/
      }
    }
  }
  //-------------------------------------------------------------------------------------------------------------------------
  
  return 0;

}

int plotPhotonPtDependence_MC(){
  
  TeresaPlottingStyle::init();


  TH1D *histoAlpha, *histoPhotonPt[nAlphaBins], *histogenJetPt[nAlphaBins];   
  TCanvas *canvas1[nPtBins][nEtaBins];
  TFile *file;
  TString fileName;
  TLatex*  info;
      
  double *x ;
  double *y;
  double *xE;
  double *yE;

  TGraphErrors* dependence[nPtBins][nEtaBins];
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Full Resolution

  for(int j=0; j<1; j++){
    for(int i=0; i<nPtBins; i++){

      x = new double[6];
      y = new double[6];
      xE = new double[6];
      yE = new double[6];


      bool first     = true;
      double norm    = 0;
      double normErr = 0;

      for(int k=0; k<nAlphaBins; k++){

	file = 0;	
	fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/hAlpha_photon_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",i+1,j+1,k+1);
	file = TFile::Open(fileName);

	
	if(file == 0){
	  cout<<"no File"<<endl;
	  continue;
	}
	
	file        -> GetObject("histo",histoAlpha);
	histoAlpha  -> SetDirectory(0); 
	
	delete file; 
	x[k]  = histoAlpha->GetMean();
	xE[k] = histoAlpha->GetMeanError();
		
	file = 0;
	//fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/hPtLeadingPhoton_intrinsic_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",i+1,j+1,k+1);
	fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/hPtLeadingJet_intrinsic_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",i+1,j+1,k+1);
	file = TFile::Open(fileName);
	
	if(file == 0){
	  cout<<"no File"<<endl;
	  continue;
	}
	 
	file             -> GetObject("histo",histoPhotonPt[k]);
	histoPhotonPt[k] -> SetDirectory(0); 
	
	

	if(first && histoPhotonPt[k]->GetMean() != 0){
	  norm    = histoPhotonPt[k]->GetMean();
	  normErr = histoPhotonPt[k]->GetMeanError();
          first = false;	  
	}

	y[k]  = (histoPhotonPt[k]->GetMean())/norm;
	yE[k] = sqrt(pow(1./norm*histoPhotonPt[k]->GetMeanError(),2)+pow(histoPhotonPt[k]->GetMean()/pow(norm,2)*normErr,2));
	
      }

     
      dependence[i][j] = new TGraphErrors(6,x,y,xE,yE);
	
      dependence[i][j] -> SetMinimum(0.8);
      dependence[i][j] -> SetMaximum(1.1);

      canvas1[i][j] = new TCanvas(Form("canvas0%i%i",i,j),"canvas",0,0,500,500);
      canvas1[i][j] -> cd();

      dependence[i][j]->Draw("AP");
	
      // Draw info boxes
      info   = new TLatex();
      info -> SetTextFont(132);
      info -> SetNDC();
      info -> SetTextSize(0.041);
	
      if(i<11) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[i],ptBins[i+1]));
      else     info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[i]));
      info->DrawLatex(0.50,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[j],etaBins[j+1]));
	
      fileName = Form("plots/PhotonDependence_intrinsic_in_%i_Pt_bin_%i_eta_bin_bin_PFCHS_mc.pdf",i,j);
      //canvas1[i][j]->SaveAs(fileName);

      delete x;
      delete y;
      delete xE;
      delete yE;
      
    }
  }
  //-------------------------------------------------------------------------------------------------------------------------
  
  return 0;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotClosure(){
  
  gStyle -> SetTitleSize(0.05,"X"); 
  gStyle -> SetTitleSize(0.05,"Y"); 
  gStyle -> SetLabelSize(0.04,"X");
  gStyle -> SetLabelSize(0.04,"Y");
  gStyle -> SetPadBottomMargin(0.15);
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.2,"X");
  gStyle -> SetTitleOffset(1.2,"Y");
  gStyle->SetMarkerStyle(20);
  
  gStyle -> SetHistLineWidth(1);
  gROOT->ForceStyle(); 
  gStyle->SetOptStat("");

  TGraphErrors *imb, *intr, *full;
  TCanvas* canvas1;
  TFile *file;
  TString fileName;
  TLatex*  info;
      
  
  int etaB = 1;
  int ptB = 3;
  
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Full Resolution

  canvas1 = new TCanvas("canvas1","canvas",0,0,500,500);
  canvas1 ->cd();

  fileName = Form("root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_imbalance_PFCHS_mc_RMS99.root",etaB,ptB);
  file = TFile::Open(fileName);
  file -> GetObject("Graph",imb);
  delete file; 

  fileName = Form("root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_intrinsic_PFCHS_mc_RMS99.root",etaB,ptB);
  file = TFile::Open(fileName);
  file        -> GetObject("Graph",intr);
  delete file; 

  fileName = Form("root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_PFCHS_mc_RMS99.root",etaB,ptB);
  file = TFile::Open(fileName);
  file        -> GetObject("Graph",full);
  delete file; 

  
  double *imbY  = imb->GetY();
  double *imbEY = imb->GetEY();

  
  double *intrY  = intr->GetY();
  double *intrEY = intr->GetEY();

  double *fullX  = full->GetX();
  double *fullEX = full->GetEX();
  double *fullY  = full->GetY();
  double *fullEY = full->GetEY();


  double *calcFullX  = new double[full->GetN()];
  double *calcFullY  = new double[full->GetN()];
  double *calcFullEX = new double[full->GetN()];
  double *calcFullEY = new double[full->GetN()];
 
 

  for(int i=0; i<full->GetN(); i++){

    calcFullX[i]  = fullX[i];
    calcFullEX[i] = fullEX[i];
    calcFullY[i]  = sqrt(imbY[i]*imbY[i] + intrY[i]*intrY[i]);
    cout<<"fullEY[i] = "<<fullEY[i]<<endl;
    cout<<"intrEY[i] = "<<intrEY[i]<<endl;
    cout<<"imbEY[i] = "<<imbEY[i]<<endl<<endl;
    calcFullY[i]  = sqrt(fullY[i]*fullY[i] - intrY[i]*intrY[i]);
    calcFullY[i]  = sqrt(imbY[i]*imbY[i] + intrY[i]*intrY[i]);
    calcFullEY[i] = 1./pow(calcFullY[i],2) *2.*(pow(imbY[i]*imbEY[i],2) + pow(intrY[i]*intrEY[i],2));
  }

  TGraphErrors *calcFull = new TGraphErrors(full->GetN(),calcFullX,calcFullY,calcFullEX,calcFullEY);
  calcFull->SetMarkerColor(1);
  calcFull->SetLineColor(1);


  //TF1* fResolutionAlpha = new TF1("fResolutionAlpha","[0] + [1]*x",0,20); 
  TF1* fResolutionAlpha = new TF1("fResolutionAlpha","TMath::Sqrt(TMath::Power([0],2) +TMath::Power([1],2) +2*[1]*[2]*x + TMath::Power(([2]*x),2) )",0,alphaBins[nAlphaBins]); 
  fResolutionAlpha -> FixParameter(1,imb->GetFunction("fResolutionAlpha")->GetParameter(0));
  fResolutionAlpha -> SetParameter(0,full->GetFunction("fResolutionAlpha")->GetParameter(0));
  fResolutionAlpha -> FixParameter(2,imb->GetFunction("fResolutionAlpha")->GetParameter(1));
  calcFull -> Fit("fResolutionAlpha","R");
  calcFull->GetFunction("fResolutionAlpha")->SetLineColor(1);

  calcFull->GetXaxis()->SetLimits(0,20);
  calcFull->SetMinimum(0.0);
  calcFull->SetMaximum(0.5);
  
  
  calcFull->Draw("AP");
  //imb->Draw("Psame");
  //full->Draw("Psame");

  

  
  full -> SetMarkerColor(2);
  full -> SetLineColor(2);
  full -> Draw("Psame");

  intr -> SetMarkerColor(3);
  intr -> SetLineColor(3);
  intr -> Draw("Psame");


  // Draw info boxes
  info   = new TLatex();
  info -> SetTextFont(132);
  info -> SetNDC();
  info -> SetTextSize(0.041);
  
  if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.62,0.80, Form(" %4.1f < #eta^{1st jet} < %4.1f",etaBins[etaB-1],etaBins[etaB]));
 
  cout<<endl;

  cout<<"real intrinsic = "<<intr->GetFunction("fResolutionAlpha")->GetParameter(0)<<endl;  
  cout<<"pred. intrinsic from full = "<<full->GetFunction("fResolutionAlpha")->GetParameter(0)<<endl;
  cout<<"pred. intrinsic from calcfull = "<<calcFull->GetFunction("fResolutionAlpha")->GetParameter(0)<<endl;
  cout<<"real imbalance = "<<imb->GetFunction("fResolutionAlpha")->GetParameter(0)<<endl;
  cout<<"imbalance from full = "<<full->GetFunction("fResolutionAlpha")->GetParameter(1)<<endl;
  cout<<"real imb**2 + real intr**2 = "<<sqrt(pow(intr->GetFunction("fResolutionAlpha")->GetParameter(0),2) + pow(imb->GetFunction("fResolutionAlpha")->GetParameter(0),2))<<endl;
  cout<<"y-intercept full = "<<sqrt(pow(full->GetFunction("fResolutionAlpha")->GetParameter(0),2) + pow(full->GetFunction("fResolutionAlpha")->GetParameter(1),2))<<endl;
  cout<<"rel error = "<<full->GetFunction("fResolutionAlpha")->GetParameter(0)/intr->GetFunction("fResolutionAlpha")->GetParameter(0) - 1.<<endl;
  cout<<"rel error bzgl. calcFull= "<<calcFull->GetFunction("fResolutionAlpha")->GetParameter(0)/intr->GetFunction("fResolutionAlpha")->GetParameter(0) - 1.<<endl;
  cout<<"chi2/ndof = "<<intr->GetFunction("fResolutionAlpha")->GetChisquare()<<"/"<<intr->GetFunction("fResolutionAlpha")->GetNDF()<<endl;



  canvas1->SaveAs("plots/fullResponseExample.pdf");




  cout<<endl<<"TEST:"<<endl;


  TF1* ftogether = new TF1("ftogether","TMath::Sqrt(TMath::Power([0],2) +TMath::Power([1],2) +2*[1]*[2]*x + TMath::Power(([2]*x),2) )",0,alphaBins[nAlphaBins]); 
  ftogether -> FixParameter(0,0.1);
  ftogether -> FixParameter(1,0.5);
  ftogether -> FixParameter(2,0.22);
  //double *togetherX =


  TF1* fhorizontal = new TF1("fhorizontal","[0]"); 
  fhorizontal ->FixParameter(0,0.1);

  return 0;
  //-------------------------------------------------------------------------------------------------------------------------
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotDataMCComparison(){
  
  TH1D *histoMC, *histoData; 
  
  TCanvas *canvas;
  TFile   *file;
  TString fileName;
  TLegend *leg_hist;
  
  gStyle->SetOptStat(0);
       
  //const int nPtBins = 12;
  //const double ptBins[pt_int+1]          = {22,36,60,88,105,148.5,165,176,200,250,300,400,1000}; 
  TString histoName = "";


  histoName = "response_photon_in_11_Pt_bin_1_eta_bin_1_alpha_bin";

  for(int i = 1; i<2; i++){

    fileName = (TString) "root_files/" + histoName + (TString) "_PFCHS_mc.root"; 
    file = new TFile(fileName);      
    histoMC =  (TH1D*) gDirectory->Get("histo");
    histoMC -> SetDirectory(0); 
    delete file; 
    
    fileName = (TString) "../data/root_files/" + histoName + (TString) "_PFCHS_data.root"; 
    file = new TFile(fileName);      
    histoData =  (TH1D*) gDirectory->Get("histo");
    histoData -> SetDirectory(0);
    delete file; 

    histoMC   -> Rebin(50);
    histoData -> Rebin(50);
    
    //histoMC   -> Sumw2();
    //histoData -> Sumw2();


    canvas = new TCanvas("canvas","canvas",100,100,500,500);
    canvas -> cd();
    //canvas -> SetLogy();
    //canvas -> SetLogx();

    //histoMC   -> GetXaxis() -> SetRange(ptBins[i]+10,ptBins[i+1]);
    //histoData -> GetXaxis() -> SetRange(ptBins[i]+10,ptBins[i+1]);

    //histoMC   -> GetXaxis()->SetRangeUser(0,0.1);  
    //histoData -> GetXaxis()->SetRangeUser(0,0.1); 

    histoMC -> Scale(histoData->Integral()/histoMC->Integral()); 
    //histoMC->SetMinimum(0);
    histoMC->SetMaximum(histoData->GetMaximum()*1.5);
    
    //histoMC -> Divide(histoMC,histoData);

    
    leg_hist = new TLegend(0.1,0.7,0.9,0.9);
    leg_hist -> SetTextSize(0.033);
    leg_hist -> SetFillColor(0);
    
    histoMC  -> SetLineColor(2);
    histoMC  -> SetMarkerColor(2);
    histoMC  -> SetMarkerStyle(20);
    TString entry = Form("MC (n = %i, #splitline{mean = %4.3f #pm %4.3f}{RMS = %4.3f #pm %4.3f}",(int)histoMC->GetEntries(),histoMC->GetMean(),histoMC->GetMeanError(),histoMC->GetRMS(),histoMC->GetRMSError());
    leg_hist -> AddEntry(histoMC,entry,"pf");

     

    histoMC->Draw("");

    histoData -> SetLineColor(3);
    histoData -> SetMarkerColor(3);
    histoData -> SetMarkerStyle(20);
    entry = Form("Data (n = %i, #splitline{mean = %4.3f #pm %4.3f}{RMS = %4.3f #pm %4.3f}",(int)histoData->GetEntries(),histoData->GetMean(),histoData->GetMeanError(),histoData->GetRMS(),histoData->GetRMSError());
    leg_hist -> AddEntry(histoData,entry,"pf");
    histoData->Draw("same");
    
    leg_hist->Draw("same");


    //canvas ->SetLogy();  
    fileName = (TString) "histo_" + histoName + (TString) ".pdf";
    canvas->SaveAs(fileName);

    //delete canvas;
    //delete leg_hist;
    //delete histoMC;
    //delete histoData;

  }

  return 0;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotNPU(){
  
  TH1D *histo, *histo1; 
  
  TCanvas* canvas;
  TFile *file;
  TString fileName;
  //TLatex*  info;
  
  //gStyle->SetOptStat("emr");
      
  canvas = new TCanvas("canvas","canvas",0,0,500,500);
  canvas ->cd();
 
  file = new TFile("root_files/hNPU.root");      
  histo =  (TH1D*) gDirectory->Get("histo");
  histo -> SetDirectory(0);
  delete file; 

  file = new TFile("root_files/hNPV.root");      
  histo1 =  (TH1D*) gDirectory->Get("histo");
  histo1 -> SetDirectory(0);
  delete file; 

  TLegend* leg_hist       =  0 ;
  leg_hist = new TLegend(0.5,0.75,0.9,0.9);
  leg_hist->SetTextSize(0.033);
  leg_hist -> SetFillColor(0);
  
  histo -> SetLineColor(2);
  histo -> SetMarkerColor(2);
  leg_hist -> AddEntry(histo,"NPU");

  histo->Draw("");

  histo1 -> SetLineColor(3);
  histo1 -> SetMarkerColor(3);
  leg_hist -> AddEntry(histo1,"NPV");
  histo1->Draw("same");

  leg_hist->Draw("same");

  canvas ->SetLogy();
  
  canvas->SaveAs("histo.pdf");

  return 0;

}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int compareSigma(){

 
  
  TGraphErrors *sigma2, *sigma3, *sigma4, *sigma5;
  
  
  TCanvas* canvas[13];
  TFile *file;
  TString fileName;

  //char form[100] = "intrinsic_";
  char form[100] = "imbalance_";
  //char form[100] = "";

 

  for(int i=0; i<1; i++){

    fileName.Form("c%i",i);
    canvas[i] = new TCanvas(fileName,fileName,0,0,500,500);
    
    canvas[i] ->cd();
 
    fileName.Form("root_files_2sigma/Resolution_for_%i_eta_bin_%sPFCHS_mc.root",i+1,form);  
    file = TFile::Open(fileName);
    file -> GetObject("Graph",sigma2);
    sigma2 -> SetMinimum(0.00);
    sigma2 -> SetMaximum(0.04);
    sigma2 -> SetMarkerColor(3);
    sigma2 -> SetLineColor(3);
    sigma2 -> GetFunction("fResolution")->SetLineColor(3);
    sigma2 -> Draw("AP");     
    delete file;

    fileName.Form("root_files_3sigma/Resolution_for_%i_eta_bin_%sPFCHS_mc.root",i+1,form);  
    file = TFile::Open(fileName);
    file -> GetObject("Graph",sigma3);
    sigma3 -> SetMarkerColor(4);
    sigma3 -> SetLineColor(4);
    sigma3 -> GetFunction("fResolution")->SetLineColor(4);
    sigma3 -> Draw("Psame");     
    delete file;

    fileName.Form("root_files_4sigma/Resolution_for_%i_eta_bin_%sPFCHS_mc.root",i+1,form);  
    file = TFile::Open(fileName);
    file -> GetObject("Graph",sigma4);
    sigma4 -> SetMarkerColor(1);
    sigma4 -> SetLineColor(1);
    sigma4 -> GetFunction("fResolution")->SetLineColor(1);
    sigma4 -> Draw("Psame");     
    delete file;

    fileName.Form("root_files_5sigma/Resolution_for_%i_eta_bin_%sPFCHS_mc.root",i+1,form);  
    file = TFile::Open(fileName);
    file -> GetObject("Graph",sigma5);
    sigma5 -> SetMarkerColor(2);
    sigma5 -> SetLineColor(2);
    sigma5 -> GetFunction("fResolution")->SetLineColor(2);
    sigma5 -> Draw("Psame");     
    delete file;

    TLegend *legend  = new TLegend(0.6,0.7,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(sigma2,"2 #sigma Range","pfl");
    legend -> AddEntry(sigma3,"3 #sigma Range","pfl");
    legend -> AddEntry(sigma4,"4 #sigma Range","pfl");
    legend -> AddEntry(sigma5,"5 #sigma Range","pfl");

    legend->Draw("same");

    fileName.Form("Sigma_comparison_in_%i_etaBin_%scomparison.pdf",i+1,form);
    canvas[i] -> SaveAs(fileName);
  }


  return 0;

}




//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int Comparison(){
  
  TGraphErrors* intr, *full, *imb;
  TGraphErrors* Ratio;
  
  
  TCanvas* canvas[4];
  TCanvas* canvas2[4];
  TFile *file;
  TString fileName;

  for(int j=0; j<1; j++){

    fileName.Form("c%i",j);
    canvas[j] = new TCanvas(fileName,fileName,0,0,500,500);

    canvas[j] ->cd();

    fileName.Form("root_files/Resolution_for_%i_eta_bin_PFCHS_mc.root",j+1);
    //fileName.Form("root_files/Scale_for_%i_eta_bin_PFCHS_mc.root",j+1);
    file = TFile::Open(fileName);
    file->GetObject("Graph",full);
    full -> SetMarkerColor(3);
    full -> SetLineColor(3);
    full -> GetFunction("fResolution")->SetLineColor(3);
    full ->SetMaximum(0.15);
    full ->SetMinimum(0.);
    full ->Draw("AP");
    delete file;


    fileName.Form("root_files/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc.root",j+1);
    //fileName.Form("root_files/Scale_for_%i_eta_bin_intrinsic_PFCHS_mc.root",j+1);
    file = TFile::Open(fileName);
    file->GetObject("Graph",intr);
    intr -> SetMarkerColor(2);
    intr -> SetLineColor(2);
    
    intr -> GetFunction("fResolution")->SetLineColor(2);
    intr ->Draw("Psame"); 

    fileName.Form("root_files/Resolution_for_%i_eta_bin_imbalance_PFCHS_mc.root",j+1);
    //fileName.Form("root_files/Scale_for_%i_eta_bin_intrinsic_PFCHS_mc.root",j+1);
    file = TFile::Open(fileName);
    file->GetObject("Graph",imb);
    imb -> SetMarkerColor(4);
    imb -> SetLineColor(4);
    
    imb -> GetFunction("fResolution")->SetLineColor(4);
    imb ->Draw("Psame"); 
    
    
    
    

    TLegend *legend  = new TLegend(0.6,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(intr,"intrinsic","l");
    legend -> AddEntry(full,"full resolution","l");

    legend->Draw("same");

    fileName.Form("Resolution_intrins_full_for_%i_eta_bin.pdf",j+1);
    canvas[j]->SaveAs(fileName);



    delete file;

    const int numEntries = 9;
    double xIntrinsic[numEntries] = {0};
    double yIntrinsic[numEntries] = {0};
    double xIntrinsicError[numEntries] = {0};
    double yIntrinsicError[numEntries] = {0};

    double xImbalance[numEntries] = {0};
    double yImbalance[numEntries] = {0};
    double xImbalanceError[numEntries] = {0};
    double yImbalanceError[numEntries] = {0};

    double xFull[numEntries] = {0};
    double yFull[numEntries] = {0};
    double xFullError[numEntries] = {0};
    double yFullError[numEntries] = {0};

    double x[numEntries] = {0};
    double y[numEntries] = {0};
    double xError[numEntries] = {0};
    double yError[numEntries] = {0};

    for(int i=0; i<numEntries; i++){
      intr   -> GetPoint(i+0,xIntrinsic[i],yIntrinsic[i]);
      imb    -> GetPoint(i+0,xImbalance[i],yImbalance[i]);
      full        -> GetPoint(i+0,xFull[i],yFull[i]);

      cout<<"xFull["<<i<<"] = "<<xFull[i]<<endl;
      cout<<"xIntrinsic["<<i<<"] = "<<xIntrinsic[i]<<endl;
      cout<<"xImbalance["<<i<<"] = "<<xImbalance[i]<<endl;

      xIntrinsicError[i] = intr -> GetErrorX(i+0);
      yIntrinsicError[i] = intr -> GetErrorY(i+0);
      xImbalanceError[i] = imb  -> GetErrorX(i+0);
      yImbalanceError[i] = imb  -> GetErrorY(i+0);
      xFullError[i]      = full      -> GetErrorX(i+0);
      yFullError[i]      = full      -> GetErrorY(i+0);
      
      y[i] = TMath::Sqrt(yIntrinsic[i]*yIntrinsic[i] + yImbalance[i]*yImbalance[i]);
      x[i] = 1./2.*(xFull[i]+xIntrinsic[i]);
      
      yError[i] =  TMath::Sqrt(TMath::Power(1./y[i]*yIntrinsic[i],2)*TMath::Power(yIntrinsicError[i],2) + TMath::Power(1./y[i]*yImbalance[i],2)*TMath::Power(yImbalanceError[i],2));
      cout<<"y["<<i<<"] = "<<y[i]*100<<endl<<endl;
      cout<<"yError["<<i<<"] = "<<yError[i]*100<<endl<<endl;

      xError[i] = 1./2.*TMath::Sqrt(TMath::Power(xFullError[i],2)+TMath::Power(xIntrinsicError[i],2));
    }

    fileName.Form("c2%i",j);
    canvas2[j] = new TCanvas(fileName,fileName,0,0,500,500);
    canvas2[j] ->cd();

    Ratio = new TGraphErrors(numEntries,x,y,xError,yError);
    Ratio -> SetTitle("Ratio: (full_response/intrinsic - 1)");
    Ratio -> SetMarkerStyle(20);
    Ratio -> SetMarkerSize(0.8);
    Ratio -> SetMinimum(0.04);
    Ratio -> SetMaximum(0.14);
   
    Ratio -> Draw("AP");

    //TF1* line = new TF1("line","0.",0,600);
    //line->SetLineColor(2);
    //line->Draw("same");

    full -> Draw("Psame");

    fileName.Form("MCClosure_for_%i_eta_bin.pdf",j+1);
    canvas2[j]->SaveAs(fileName);


  }


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


int makeOtherGraph(){


  TFile *file, *file1;  
  /*
  file = TFile::Open("Photon1Pt_PF_data.root");
  file->GetObject("4",hPhotonPtData);
    
  file1 = TFile::Open("Photon1Pt_PF_mc.root");
  file1->GetObject("4",hPhotonPtMC);

  hPhotonPtMC->Scale(hPhotonPtData->Integral()/hPhotonPtMC->Integral()); 

  //TCanvas* canvas = DrawComparison(hPhotonPtData,hPhotonPtMC,"Photon Pt spectrum (data and MC) New","Lumititle","Photon Pt", 1);
  
  //canvas ->Print("PhotonPtComparison.pdf");

  
  delete file;
  delete file1; 
  */ 

  TH1D* hNVtxData;
  TH1D* hNVtxMC;

  file = TFile::Open("VtxN_PF_data.root");
  file->GetObject("hVtxN",hNVtxData);

  
  file1 = TFile::Open("VtxN_PF_mc.root");
  file1->GetObject("hVtxN",hNVtxMC);

  hNVtxMC->Scale(hNVtxData->Integral()/hNVtxMC->Integral()); 

  //canvas = DrawComparison(hNVtxData, hNVtxMC,"Number of Vertices (data and MC) New PU Dist","Lumititle","#Vtx", 1);
  //canvas ->Print("NVtxComparison.pdf");

  

  char title[100];
  TH1D* hNVtxD[8];
  TH1D* hNVtxM[8];
  TCanvas* canv[8];
  TCanvas* canvas;
  for(int i=0; i<8;i++){

  delete file;
  delete file1;
 
  sprintf(title,"VtxN%i_PF_data.root",i);
  file = TFile::Open(title);
  
  sprintf(title,"hVtxPtBinned%i",i);
  file->GetObject(title,hNVtxD[i]);
  
  sprintf(title,"VtxN%i_PF_mc.root",i);
  file1 = TFile::Open(title);

  sprintf(title,"hVtxPtBinned%i",i);
  file1->GetObject(title,hNVtxM[i]);
  
  
  hNVtxM[i]->Scale(hNVtxD[i]->Integral()/hNVtxM[i]->Integral()); 
  
  sprintf(title,"Number of Vertices (data and MC) pt %i  NEW PU Dist",i);
  canv[i] = DrawComparison(hNVtxD[i], hNVtxM[i],title,"","#Vtx", 1);
  
  sprintf(title,"NVtxComparison%i.pdf",i);
  canv[i] ->Print(title);

  }

  delete file;
  delete file1;
  TH1D* hIsoData;
  TH1D* hIsoMC;
 
  sprintf(title,"hPhotonIsoEcal_PF_data.root");
  file = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoEcal");
  file->GetObject(title,hIsoData);
  
  sprintf(title,"hPhotonIsoEcal_PF_mc.root");
  file1 = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoEcal");
  file1->GetObject(title,hIsoMC);
  
  
  hIsoMC->Scale(hIsoData->Integral()/hIsoMC->Integral()); 
  
  sprintf(title,"Ecal Isolation");
  canvas = DrawComparison(hIsoData,hIsoMC,title,"","Iso", 1);
  
  sprintf(title,"PhotonIsoEcal.pdf");
  canvas ->Print(title);

  delete file;
  delete file1;
  
 
  sprintf(title,"hPhotonIsoHcal_PF_data.root");
  file = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoHcal");
  file->GetObject(title,hIsoData);
  
  sprintf(title,"hPhotonIsoHcal_PF_mc.root");
  file1 = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoHcal");
  file1->GetObject(title,hIsoMC);
  
  
  hIsoMC->Scale(hIsoData->Integral()/hIsoMC->Integral()); 
  
  sprintf(title,"Hcal Isolation");
  canvas = DrawComparison(hIsoData, hIsoMC,title,"","Iso", 1);
  
  sprintf(title,"PhotonIsoHcal.pdf");
  canvas ->Print(title);

  delete file;
  delete file1;
  
  sprintf(title,"hPhotonIsoTrk_PF_data.root");
  file = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoTrk");
  file->GetObject(title,hIsoData);
  
  sprintf(title,"hPhotonIsoTrk_PF_mc.root");
  file1 = TFile::Open(title);
  
  sprintf(title,"hPhotonIsoTrk");
  file1->GetObject(title,hIsoMC);
  
  
  hIsoMC->Scale(hIsoData->Integral()/hIsoMC->Integral()); 
  
  sprintf(title,"Trk Isolation");
  canvas = DrawComparison(hIsoData, hIsoMC,title,"","Iso", 1);
  
  sprintf(title,"PhotonIsoTrk.pdf");
  canvas ->Print(title);
  

  


  return 0;
}

int makeRhoGraph(){

  char title[100];
  TH1D* hRhoD[8];
  TH1D* hRhoM[8];
  TCanvas* canv[8];
  TFile *file;
  TFile *file1;
  for(int i=0; i<8;i++){

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
  TH1D* hTriggerBefore[8];
  TH1D* hTriggerAfter[8];
  TH1D* hTrigger[8];
  TCanvas* canv[8];
  TFile *file;
  TFile *file1;
  for(int i=0; i<8;i++){


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
   double YRangeMax = MaxY;
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
   //pad1->SetLogy();
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
   r->Draw("ep");
   TLine l;
   l.DrawLine(MinX, 0., MaxX, 0.);
   c->cd();
   return c;
}
