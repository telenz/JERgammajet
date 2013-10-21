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
#include "TLine.h"
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
#include "TROOT.h"
#include "../CODE/myDeclarations.h"
#include "plotStyle.h"
#include "/afs/naf.desy.de/user/t/telenz/comparison/tdrstyle_mod.C"

TGraphErrors* GetTGraphErrors(TString filename, TString objectName);
TF1* GetTF1(TString filename, TString objectName);
TH1D* GetTH1D(TString filename, TString objectName);
TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool isData);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// MC Closure

int MCClosure(){

  
  TeresaPlottingStyle::init();
  gStyle -> SetPadLeftMargin(0.20);
  gStyle -> SetTitleOffset(1.7,"Y");

  TString method = "RMS99";
 
  
  TGraphErrors* intr, *full, *Ratio;  
  TCanvas *canvas[nEtaBins], *canvas2[nEtaBins];
  TFile *file;
  TString fileName;
  TString etaRegion; 



  for(int j=0; j<nEtaBins; j++){

    if(j == 0) etaRegion.Form("|#eta^{1st jet}| < %4.1f", etaBins[j+1]);
    else       etaRegion.Form("%4.1f <|#eta^{1st jet}|< %4.1f",etaBins[j], etaBins[j+1]);

    fileName.Form("c%i",j);
    canvas[j] = new TCanvas(fileName,fileName,0,0,500,500);

    canvas[j] ->cd();
  
    fileName = "../plots_2012/PF_L1CHS/mc/root_files" + (TString) "/Resolution_for_" + (long) (j+1) + "_eta_bin_PFCHS_mc_" + method + ".root";    
    file =  TFile::Open(fileName);
    file -> GetObject("Graph",full);
    delete file;

    fileName = "../plots_2012/PF_L1CHS/mc/root_files" + (TString) "/Resolution_for_" + (long) (j+1) + "_eta_bin_intrinsic_PFCHS_mc_" + method + ".root"; 
    file = TFile::Open(fileName);
    file->GetObject("Graph",intr);
    delete file;

    full -> SetMarkerColor(3);
    full -> SetLineColor(3);
    full -> GetFunction("fResolution")->SetLineColor(3);
    intr -> SetMarkerColor(2);
    intr -> SetLineColor(2);
    intr -> GetFunction("fResolution")->SetLineColor(2);

    
    full ->SetMaximum(0.12);
    full ->SetMinimum(0.04);
    full ->Draw("AP");
         
    intr ->Draw("Psame"); 
    
 
    TLegend *legend  = new TLegend(0.6,0.7,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextFont(132);
    legend -> SetTextSize(0.042);
    legend -> AddEntry(intr,"intrinsic","l");
    legend -> AddEntry(full,"full resolution","l");

    legend->Draw("same");

    fileName = (TString) "plots/Resolution_intrinsic_full_for_"  + (long) (j+1) + "_eta_bin" + "_" + method + ".pdf";
    canvas[j]->SaveAs(fileName);
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    const int numEntries = intr ->GetN();
    double *xIntrinsic   = intr -> GetX();
    double *yIntrinsic   = intr -> GetY();
    double *xIntrinsicError = intr -> GetEX();
    double *yIntrinsicError  = intr -> GetEY();
    double *xFull = full -> GetX();
    double *yFull = full -> GetY();
    double *xFullError = full -> GetEX();
    double *yFullError = full -> GetEY();

    double *x = new double[numEntries];
    double *y = new double[numEntries];
    double *xError = new double[numEntries];
    double *yError = new double[numEntries];

    for(int i=0; i<numEntries; i++){
      
      
      y[i] = yFull[i]/yIntrinsic[i] -1 ;
      x[i] = 1./2.*(xFull[i]+xIntrinsic[i]);
      
      yError[i] =  TMath::Sqrt(TMath::Power((1./yIntrinsic[i]),2)*TMath::Power(yFullError[i],2)+TMath::Power((yFull[i]/(TMath::Power(yIntrinsic[i],2))),2)*TMath::Power(yIntrinsicError[i],2)-2.0*yFull[i]/TMath::Power(yIntrinsic[i],3)*yFullError[i]*yIntrinsicError[i]);
      //cout<<"y["<<i<<"] = "<<y[i]*100<<endl<<endl;
      //cout<<"yError["<<i<<"] = "<<yError[i]*100<<endl<<endl;

      xError[i] = 1./2.*TMath::Sqrt(TMath::Power(xFullError[i],2)+TMath::Power(xIntrinsicError[i],2));
  }

    fileName.Form("c2%i",j);
    canvas2[j] = new TCanvas(fileName,fileName,0,0,500,500);
    canvas2[j] ->cd();

    Ratio = new TGraphErrors(numEntries,x,y,xError,yError);

       
    Ratio -> SetTitle("");
    Ratio -> GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    Ratio -> GetXaxis()->SetLimits(0,600);
    Ratio -> GetYaxis()->SetTitle("#frac{predic. intrinsic}{intrinsic} -1");
    Ratio -> SetMaximum(0.3);
    Ratio -> SetMinimum(-0.3);
    Ratio -> Draw("AP");

    TF1* line = new TF1("line","0.",0,600);
    
    line->SetLineColor(2);
    line->Draw("same");

    TLatex*  info   = new TLatex();
    info-> SetNDC();
    info->SetTextSize(0.05);
    info->DrawLatex(0.60, 0.80,etaRegion);
    

    fileName =  (TString) "plots/MCClosure_for_" + (long) (j+1) + "_eta_bin" + "_" + method + ".pdf";
    canvas2[j]->SaveAs(fileName);
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
  }

  return 0;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotJetResponseAlpha(){


  gDirectory->Delete(); 
  gROOT->GetListOfCanvases()->Delete();
  
  TeresaPlottingStyle::init();


  bool addAll = true;
  bool addMC  = true;



  TString sourceResponseAlpha, sourceIntrinsicAlpha, sourceImbalanceAlpha, sourceTotalAlpha, sourceData, pdfFile, etaRegion, ptRegion, title, EtaPtRegion; 
    
  TLatex*  info;
  TLegend *legend;
  TGraphErrors* Add;
  TF1 *fScaleAlpha;
  TF1 *total;
  for(int j=0; j<nEtaBins; j++){
    for(int i=0;i<nPtBins;i++){       


      sourceResponseAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_PFCHS_mc_RMS99.root",j+1,i+1);
      sourceIntrinsicAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_intrinsic_PFCHS_mc_RMS99.root",j+1,i+1);
      sourceImbalanceAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_imbalance_PFCHS_mc_RMS99.root",j+1,i+1);
      sourceTotalAlpha.Form("../plots_2012/PF_L1CHS/mc/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_total_PFCHS_mc_RMS99.root",j+1,i+1);
      sourceData.Form("../plots_2012/PF_L1CHS/data/root_files/jet_energy_resolution_for_%i_eta_bin_%i_pTGamma_bin_PFCHS_data_RMS99.root",j+1,i+1);
      
      if(addAll)     pdfFile.Form("plots/JER_for_%i_eta_bin_%i_pTGamma_bin_all_contributions_PFCHS_RMS99_mc.pdf",j+1,i+1);
      else if(addMC) pdfFile.Form("plots/JER_for_%i_eta_bin_%i_pTGamma_bin_wo_data_PFCHS_RMS99_mc.pdf",j+1,i+1);
      else           pdfFile.Form("plots/JER_for_%i_eta_bin_%i_pTGamma_bin_wo_data_mc_PFCHS_RMS99_mc.pdf",j+1,i+1);
      title = "Jet Energy Resolution";
 
     

      TCanvas *c = new TCanvas(EtaPtRegion,title,200,10,450,450);
      c -> cd();
      
      legend  = new TLegend(0.45,0.65,0.9,0.9);
      legend -> SetFillColor(0);
      legend -> SetTextFont(132);
      legend -> SetTextSize(0.042);
  
      
  
      // Intrinsic
      Add = GetTGraphErrors(sourceIntrinsicAlpha,"Graph");  
      if(Add == 0) continue;
      Add -> SetMinimum(0.0);
      Add -> SetMaximum(0.4);   
      Add -> GetXaxis()->SetLimits(0,20);
      fScaleAlpha = Add->GetFunction("fResolutionAlpha");
  
      Add -> SetMarkerColor(46);
      Add -> SetLineColor(46);
      fScaleAlpha -> SetLineColor(46);
      Add -> GetXaxis() -> SetTitle("#alpha [%]");   
      Add -> GetXaxis() -> SetTitle("p_{T}^{2nd jet}/p_{T}^{#gamma} #upoint 100");    
      Add -> GetYaxis() -> SetTitle("Resolution"); 
      Add -> SetTitle("");
      legend -> AddEntry(Add,"Intrinsic","l");
      Add -> Draw("AP");

      Add = GetTGraphErrors(sourceImbalanceAlpha,"Graph");  
      fScaleAlpha = Add->GetFunction("fResolutionAlpha");
  
      Add -> SetMarkerColor(4);
      Add -> SetLineColor(4);
      fScaleAlpha -> SetLineColor(4);
      legend -> AddEntry(Add,"Imbalance","l");
      Add -> Draw("Psame");


      total = GetTF1(sourceTotalAlpha,"totalResolution");
      total -> SetLineColor(14);
      total -> SetLineWidth(3);
      legend -> AddEntry(total,"Total","l");
      total -> Draw("same");


      // MC pseudo data
      Add = GetTGraphErrors(sourceResponseAlpha,"Graph");  
      fScaleAlpha = Add->GetFunction("fResolutionAlpha");
      Add -> SetMarkerColor(1);
      Add -> SetLineColor(1);
      fScaleAlpha -> SetLineColor(1);
      Add -> GetXaxis() -> SetTitle("#alpha [%]");   
      Add -> GetXaxis() -> SetTitle("p_{T}^{2nd jet}/p_{T}^{#gamma} #upoint 100");    
      Add -> GetYaxis() -> SetTitle("Resolution"); 
      Add -> SetTitle("");           
      if(addMC){
	legend -> AddEntry(Add,"#gamma + Jet (pseudo data)","l");
	Add    -> Draw("Psame");
	total  -> Draw("same");
      }

      // Data
      Add = GetTGraphErrors(sourceData,"Graph");  
      if(Add != 0){
	fScaleAlpha = Add->GetFunction("fResolutionAlpha");
	Add -> SetMarkerColor(8);
	Add -> SetLineColor(8);
	fScaleAlpha -> SetLineColor(8);
	if(addAll){
	  legend -> AddEntry(Add,"Data","l");
	  Add -> Draw("Psame");
	}
      }

      // Draw info boxes
      if (i==nPtBins-1) ptRegion.Form(" %4.1f GeV < p_{T}^{#gamma} \n",ptBins[i]);
      else       ptRegion.Form("%4.1f GeV < p_{T}^{#gamma} < %4.1f GeV \n" ,ptBins[i], ptBins[i+1]);
      
      if(j == 0) etaRegion.Form("|#eta^{1st jet}| < %4.1f", etaBins[j+1]);
      else       etaRegion.Form("%4.1f < |#eta^{1st jet}| < %4.1f",etaBins[j], etaBins[j+1]);
      
      EtaPtRegion.Form("#splitline{%s}{%s}",ptRegion.Data(),etaRegion.Data());

      info   = new TLatex();
  
      info -> SetNDC();    
      info->DrawLatex(0.18,0.60,  ptRegion);
      info->DrawLatex(0.18,0.53,  etaRegion);
      info->DrawLatex(0.50,0.53, "Anti-k_{T} 0.5 PFCHSJets");
      
      
      legend -> Draw("same");
      c -> SaveAs(pdfFile);
    }
  }

 
  return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotResponse(){
  
 
  TeresaPlottingStyle::init();

  gStyle -> SetOptLogy();
  gStyle -> SetMarkerStyle(1);



  TH1D *histoPhoton1, *histoJet1, *histoPhoton2 , *histoJet2, *histoJet3, *histoTogether1, *histoTogether2;   
  TCanvas* canvas1, *canvas2, *canvas3, *canvas4, *canvas5;
  TString fileName;
  TLatex*  info;
      
  
  const int etaB = 1;
  const int alphaB = 6;
  const int ptB = 12;
  
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Full Resolution

  canvas1 = new TCanvas("canvas1","canvas",0,0,500,500);
  canvas1 ->cd();

  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_photon_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoPhoton1 = GetTH1D(fileName,"histo");
  
  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoJet1 = GetTH1D(fileName,"histo");

  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoTogether1 = GetTH1D(fileName,"histo");

  histoTogether1 -> Add(histoPhoton1);
  histoTogether1 -> Scale(histoTogether1->GetEntries()/histoTogether1->Integral());
  histoTogether1 -> Rebin(40);
  histoTogether1 -> GetXaxis() -> SetRangeUser(0.,2.);

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
  info -> SetNDC();
  
  if(ptB<12) info->DrawLatex(0.58,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.58,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.58,0.78, Form(" %4.1f < |#eta^{1st jet}| < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
  
  canvas1->SaveAs("plots/fullResponseExample.pdf");
  //-------------------------------------------------------------------------------------------------------------------------

  canvas4 = new TCanvas("canvas4","canvas",0,600,500,500);
  canvas4 ->cd();

  histoTogether1 -> Draw(); 
  
  double scaleFactor = histoJet1->Integral() + histoPhoton1->Integral(); 

  histoJet1 -> Scale(histoTogether1->GetEntries()/scaleFactor);

  histoJet1 -> Rebin(40);
  histoJet1 -> GetXaxis() -> SetRangeUser(0.,2.);
  histoJet1 -> SetMarkerColor(8);
  histoJet1 -> SetLineColor(8);
  histoJet1 -> Draw("same");

  histoPhoton1 -> Scale(histoTogether1->GetEntries()/scaleFactor);
  histoPhoton1 -> Rebin(40);
  histoPhoton1 -> GetXaxis() -> SetRangeUser(0.,2.);
  histoPhoton1 -> SetMarkerColor(2);
  histoPhoton1 -> SetLineColor(2);

  histoPhoton1 -> Draw("same");


  // Draw info boxes
  info   = new TLatex();
  info -> SetNDC();

  
  //if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  //else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  //info->DrawLatex(0.62,0.80, Form(" %4.1f < |#eta^{1st jet}| < %4.1f",etaBins[etaB-1],etaBins[etaB]));
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
  
  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_imbalance_photon_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoPhoton2 = GetTH1D(fileName,"histo");
  

  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_imbalance_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoJet2 = GetTH1D(fileName,"histo");

  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_imbalance_jet_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoTogether2 = GetTH1D(fileName,"histo");

  histoTogether2 -> Add(histoPhoton2);
  histoTogether2 -> Scale(histoTogether2->GetEntries()/histoTogether2->Integral());
  histoTogether2 -> Rebin(40);
  histoTogether2 -> GetXaxis() -> SetRangeUser(0.,2.);

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
  info -> SetNDC();
  
  if(ptB<12) info->DrawLatex(0.58,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.58,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.58,0.78, Form(" %4.1f < |#eta^{1st jet}| < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
  
  canvas2->SaveAs("plots/imbalanceExample.pdf");

  //-------------------------------------------------------------------------------------------------------------------------

  canvas5 = new TCanvas("canvas5","canvas",500,600,500,500);
  canvas5 ->cd();

  histoTogether2 -> Draw(); 
  
  scaleFactor = histoJet2->Integral() + histoPhoton2->Integral(); 

  histoJet2 -> Scale(histoTogether2->GetEntries()/scaleFactor);

  histoJet2 -> Rebin(40);
  histoJet2 -> GetXaxis() -> SetRangeUser(0.,2.);
  histoJet2 -> SetMarkerColor(8);
  histoJet2 -> SetLineColor(8);
  histoJet2 -> Draw("same");

  histoPhoton2 -> Scale(histoTogether2->GetEntries()/scaleFactor);
  histoPhoton2 -> Rebin(40);
  histoPhoton2 -> GetXaxis() -> SetRangeUser(0.,2.);
  histoPhoton2 -> SetMarkerColor(2);
  histoPhoton2 -> SetLineColor(2);
  histoPhoton2 -> Draw("same");


  // Draw info boxes
  info   = new TLatex();
  info -> SetNDC();
  
  //if(ptB<12) info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  //else       info->DrawLatex(0.62,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  //info->DrawLatex(0.62,0.80, Form(" %4.1f < |#eta^{1st jet}| < %4.1f",etaBins[etaB-1],etaBins[etaB]));
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
  
  fileName = Form("../plots_2012/PF_L1CHS/mc/root_files/response_intrinsic_in_%i_Pt_bin_%i_eta_bin_%i_alpha_bin_PFCHS_mc.root",ptB,etaB,alphaB);
  histoJet3 = GetTH1D(fileName,"histo");

  histoJet3 ->Scale(histoJet3->GetEntries()/histoJet3->Integral());
  histoJet3 -> Rebin(40);
  histoJet3 -> GetXaxis() -> SetRangeUser(0.,2.);

  histoJet3 -> GetXaxis()->SetTitle("p_{T}^{reco. jet}/p_{T}^{gen. jet}");
  histoJet3 -> GetYaxis()->SetTitle("# Events");
  histoJet3 -> SetTitle("Intrinsic");
  histoJet3 -> SetMinimum(1.);
  histoJet3 -> SetMaximum(histoJet3->GetMaximum()*20.);

  histoJet3 -> SetMarkerColor(46);
  histoJet3 -> SetLineColor(46);
  
  histoJet3 -> Draw(); 


  // Draw info boxes
  info   = new TLatex();
  info -> SetNDC();
  
  if(ptB<12) info->DrawLatex(0.58,0.85, Form("%4.0f GeV< p_{T}^{#gamma} < %4.0f GeV",ptBins[ptB-1],ptBins[ptB]));
  else       info->DrawLatex(0.58,0.85, Form("%4.0f GeV< p_{T}^{#gamma} ",ptBins[ptB-1]));
  info->DrawLatex(0.58,0.78, Form(" %4.1f < |#eta^{1st jet}| < %4.1f",etaBins[etaB-1],etaBins[etaB]));
  info->DrawLatex(0.18,0.85, Form("%4.3f < #alpha < %4.3f",alphaBins[alphaB-1]/100.,alphaBins[alphaB]/100.));
  
  canvas3->SaveAs("plots/intrinsicExample.pdf");
  
  return 0;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int plotResolutionOfPhotonPt(){

  TeresaPlottingStyle::init();

  TString sourceData, sourceIntrinsic, sourceMC, sourceImbalance, pdfFile, etaRegion, ptRegion, title, EtaPtRegion; 
    
  TLatex*  info;
  TGraphErrors* intrinsicJER, *imbalanceJER, *dataJER, *mcJER;
  TCanvas *c[nEtaBins], *cImb[nEtaBins];
  for(int j=0; j<nEtaBins; j++){

    sourceIntrinsic.Form("../plots_2012/PF_L1CHS/mc/root_files/Resolution_for_%i_eta_bin_intrinsic_PFCHS_mc_RMS99.root",j+1);
    sourceImbalance.Form("../plots_2012/PF_L1CHS/mc/root_files/Resolution_for_%i_eta_bin_imbalance_PFCHS_mc_RMS99.root",j+1);
    sourceMC.Form("../plots_2012/PF_L1CHS/mc/root_files/Resolution_for_%i_eta_bin_PFCHS_mc_RMS99.root",j+1);
    sourceData.Form("../plots_2012/PF_L1CHS/data/root_files/Resolution_for_%i_eta_bin_PFCHS_data_RMS99.root",j+1);
    pdfFile.Form("plots/Resolution_for_%i_eta_bin_PFCHS_mc_RMS99.pdf",j+1);


    c[j] = new TCanvas(EtaPtRegion,title,200,10,450,450);
    TString name = "imbalance";
    cImb[j] = new TCanvas(name,name,200,10,450,450);
    c[j] -> cd();
    TLegend *legend = new TLegend(0.4,0.75,0.9,0.9);      

    intrinsicJER = GetTGraphErrors(sourceIntrinsic,"Graph");  
    imbalanceJER = GetTGraphErrors(sourceImbalance,"Graph");
    dataJER      = GetTGraphErrors(sourceData,"Graph");  
    mcJER        = GetTGraphErrors(sourceMC,"Graph");  
    if(intrinsicJER  == 0) continue;
      
    intrinsicJER -> GetFunction("fResolution")-> SetLineColor(46);
    intrinsicJER -> SetMarkerColor(46);
    intrinsicJER -> SetLineColor(46);
    imbalanceJER -> GetFunction("fResolution")-> SetBit(TF1::kNotDraw);
    imbalanceJER -> SetMarkerColor(4);
    imbalanceJER -> SetLineColor(4);
    mcJER -> GetFunction("fResolution")-> SetLineColor(1);
    mcJER -> SetMarkerColor(1);
    mcJER -> SetLineColor(1);
    dataJER -> GetFunction("fResolution")-> SetLineColor(8);
    dataJER -> SetMarkerColor(8);
    dataJER -> SetLineColor(8);


    intrinsicJER -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");   
    imbalanceJER -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
    mcJER        -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");   
    dataJER      -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");   
    intrinsicJER -> GetYaxis() -> SetTitle("Resolution"); 
    imbalanceJER -> GetYaxis() -> SetTitle("Residual Imbalance"); 
    mcJER        -> GetYaxis() -> SetTitle("Resolution"); 
    dataJER      -> GetYaxis() -> SetTitle("Resolution"); 
    intrinsicJER -> SetTitle("");
    imbalanceJER -> SetTitle("");
    mcJER        -> SetTitle("");
    dataJER      -> SetTitle("");
                        
    intrinsicJER -> SetMinimum(0.00);
    intrinsicJER -> SetMaximum(0.15);   
    intrinsicJER -> GetXaxis() -> SetRangeUser(0,600);
    imbalanceJER -> SetMinimum(0.00);
    imbalanceJER -> SetMaximum(0.15);   
    imbalanceJER -> GetXaxis() -> SetRangeUser(0,600);

    mcJER -> SetMinimum(0.00);
    mcJER -> SetMaximum(0.15);   
    mcJER -> GetXaxis() -> SetRangeUser(0,600);
      
    //intrinsicJER -> Draw("AP");
    mcJER        -> Draw("AP");
    dataJER      -> Draw("Psame");

    //legend -> SetTextSize(0.033);
    legend -> SetFillColor(0);
    legend -> SetLineWidth(2);
    legend -> AddEntry(dataJER, "Data");
    legend -> AddEntry(mcJER, "Simulation");
    legend -> Draw("same");

    // Draw info boxes
      
    if(j == 0) etaRegion.Form("|#eta^{1st jet}| < %4.1f", etaBins[j+1]);
    else       etaRegion.Form("%4.1f < |#eta^{1st jet}| < %4.1f",etaBins[j], etaBins[j+1]);
      
    info   = new TLatex();
  
    info -> SetNDC();    
    info->DrawLatex(0.53,0.65,  etaRegion);
    info->DrawLatex(0.40,0.25, "Anti-k_{T} 0.5 PFCHSJets");
      
    c[j] -> SaveAs(pdfFile);

    // Draw imbalance
    cImb[j] -> cd(); 
    imbalanceJER -> Draw("AP");
    delete info;
    info   = new TLatex();
    info -> SetNDC();    
    info->DrawLatex(0.45,0.65,  etaRegion);
    info->DrawLatex(0.45,0.75, "Anti-k_{T} 0.5 PFCHSJets");

    pdfFile.Form("plots/Imbalance_for_%i_eta_bin_PFCHS_mc_RMS99.pdf",j+1);
    cImb[j] -> SaveAs(pdfFile);
    
  }

  return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Ratio Closure

void MCClosureRatioALLMETHODS(){

  //setTDRStyle(true);
  TeresaPlottingStyle::init();
  gStyle->SetErrorX(0);
  //gStyle->SetHistLineWidth(2); 
  gStyle -> SetNdivisions(505, "X");
  gStyle -> SetNdivisions(508, "Y");

   double eta_bins[5] = {0., 0.5, 1.1, 1.7, 2.3};
   
   TH1F *Res_INPUT = new TH1F("INPUT","", 4, eta_bins);
   TH1F *Res_Gaus  = new TH1F("MCsmeared_MC_ratio_Gaus", "", 4, eta_bins);
   TH1F *Res_RMS99 = new TH1F("MCsmeared_MC_ratio_RMS99", "", 4, eta_bins);
   TH1F *Res_RMS95 = new TH1F("MCsmeared_MC_ratio_RMS95", "", 4, eta_bins);

   // Input Parameters
   Res_INPUT->SetBinContent(1, 1.05);
   Res_INPUT->SetBinContent(2, 1.07);
   Res_INPUT->SetBinContent(3, 1.09);
   Res_INPUT->SetBinContent(4, 1.11);

   Res_INPUT->SetBinError(1, 0);
   Res_INPUT->SetBinError(2, 0);
   Res_INPUT->SetBinError(3, 0);
   Res_INPUT->SetBinError(4, 0);
   
   // 2sigma Gauss CHS
   Res_Gaus->SetBinContent(1, 1.044);
   Res_Gaus->SetBinContent(2, 1.066);
   Res_Gaus->SetBinContent(3, 1.082);
   Res_Gaus->SetBinContent(4, 1.099);
   
   Res_Gaus->SetBinError(1, 0.005);
   Res_Gaus->SetBinError(2, 0.005);
   Res_Gaus->SetBinError(3, 0.007);
   Res_Gaus->SetBinError(4, 0.017);

   Res_Gaus->SetBinContent(1, 1.043);
   Res_Gaus->SetBinContent(2, 1.062);
   Res_Gaus->SetBinContent(3, 1.084);
   Res_Gaus->SetBinContent(4, 1.105);
   
   Res_Gaus->SetBinError(1, 0.005);
   Res_Gaus->SetBinError(2, 0.005);
   Res_Gaus->SetBinError(3, 0.008);
   Res_Gaus->SetBinError(4, 0.019);

   // RMS 99
   Res_RMS99->SetBinContent(1, 1.048);
   Res_RMS99->SetBinContent(2, 1.068);
   Res_RMS99->SetBinContent(3, 1.086);
   Res_RMS99->SetBinContent(4, 1.116);
   
   Res_RMS99->SetBinError(1, 0.004);
   Res_RMS99->SetBinError(2, 0.004);
   Res_RMS99->SetBinError(3, 0.005);
   Res_RMS99->SetBinError(4, 0.013);

   Res_RMS99->SetBinContent(1, 1.048);
   Res_RMS99->SetBinContent(2, 1.067);
   Res_RMS99->SetBinContent(3, 1.086);
   Res_RMS99->SetBinContent(4, 1.116);
   
   Res_RMS99->SetBinError(1, 0.004);
   Res_RMS99->SetBinError(2, 0.004);
   Res_RMS99->SetBinError(3, 0.006);
   Res_RMS99->SetBinError(4, 0.013);

   // RMS 95   
   Res_RMS95->SetBinContent(1, 1.048);
   Res_RMS95->SetBinContent(2, 1.067);
   Res_RMS95->SetBinContent(3, 1.085);
   Res_RMS95->SetBinContent(4, 1.101);
   
   Res_RMS95->SetBinError(1, 0.004);
   Res_RMS95->SetBinError(2, 0.004);
   Res_RMS95->SetBinError(3, 0.006);
   Res_RMS95->SetBinError(4, 0.014);

   Res_RMS95->SetBinContent(1, 1.049);
   Res_RMS95->SetBinContent(2, 1.067);
   Res_RMS95->SetBinContent(3, 1.085);
   Res_RMS95->SetBinContent(4, 1.098);
   
   Res_RMS95->SetBinError(1, 0.004);
   Res_RMS95->SetBinError(2, 0.004);
   Res_RMS95->SetBinError(3, 0.006);
   Res_RMS95->SetBinError(4, 0.014);
   

   ///////////////////////////////////////////////////////////////
  
   TCanvas *c = new TCanvas();
   Res_INPUT->SetXTitle("|#eta^{1st jet}|");
   Res_INPUT->GetXaxis()->SetRangeUser(0., 5.);
   Res_INPUT->SetYTitle("MC_{smeared}/MC");
   Res_INPUT->GetYaxis()->SetRangeUser(1.04, 1.14);
   Res_INPUT->SetMarkerStyle(20);
   Res_INPUT->SetMarkerSize(1.2);
   Res_INPUT->SetLineColor(kBlack);
   Res_INPUT->SetMarkerColor(kBlack);
   Res_INPUT->Draw("e1p");

   Res_Gaus->SetMarkerStyle(23);
   Res_Gaus->SetMarkerSize(1.2);
   Res_Gaus->SetLineColor(kBlue-4);
   Res_Gaus->SetMarkerColor(kBlue-4);
   Res_Gaus->DrawClone("e1psame");

   Res_INPUT->Draw("e1psame");

   Res_RMS99->SetMarkerStyle(21);
   Res_RMS99->SetMarkerSize(1.2);
   Res_RMS99->SetLineColor(kPink-8);
   Res_RMS99->SetMarkerColor(kPink-8);
   Res_RMS99->Draw("e1psame");

   Res_RMS95->SetMarkerStyle(24);
   Res_RMS95->SetMarkerSize(1.2);
   Res_RMS95->SetLineColor(kGreen+2);
   Res_RMS95->SetMarkerColor(kGreen+2);
   Res_RMS95->Draw("e1psame");
   Res_RMS95->SetMarkerSize(1.0);
   Res_RMS95->Draw("e1psame");

   //cmsPrel();

   TLegend *leg = new TLegend(0.20, 0.70, 0.40, 0.90);
   leg->SetBorderSize(0);
   // leg->SetBorderMode(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.045);

   leg->AddEntry(Res_INPUT,"Input numbers", "pfl");
   leg->AddEntry(Res_Gaus, "Gaus fit 2 #sigma interval", "pfl");
   leg->AddEntry(Res_RMS99,"St. dev. 99%", "pfl");
   leg->AddEntry(Res_RMS95,"St. dev. 95%", "pfl");

   leg->Draw("same");
   c->Print("plots/MCClosureRatio.pdf");   
   //c->Print("MCClosure.pdf");   
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int drawDeltaPhiCut(bool isData = true){

  TeresaPlottingStyle::init();

  TFile *file;
  TString filename;
  TString saveName;
  TH1D* histo;  
  TCanvas* c = new TCanvas("c","c",200,10,450,450);
  c -> cd();
  c -> SetLogy(); 
 
  if(isData) filename = "../plots_2012/PF_L1CHS/data/root_files/hDeltaPhi_PFCHS_data.root";
  else       filename = "../plots_2012/PF_L1CHS/mc/root_files/hDeltaPhi_PFCHS_mc.root";

  file = TFile::Open(filename);  
  file -> GetObject("histo",histo);
  histo-> SetTitle("#Delta #Phi after all cuts");
  histo -> GetYaxis()->SetRangeUser(1,histo->GetMaximum()*5);
  histo -> GetXaxis()->SetRangeUser(1.5,3.2);
  histo -> Draw();

  TLine *line = new TLine(2.95,histo->GetMinimum(),2.95,histo->GetMaximum());
  line -> SetLineColor(2);
  line -> SetLineWidth(2);
  line -> Draw("same");



  TLatex *info   = new TLatex();
  
  info -> SetNDC();    
  info -> SetTextSize(0.06);
  if(isData) info->DrawLatex(0.20,0.60, "Data");
  else       info->DrawLatex(0.20,0.60, "MC");

  if(isData) saveName = "plots/hDeltaPhi_PFCHS_data.pdf";
  else       saveName = "plots/hDeltaPhi_PFCHS_mc.pdf";

  c -> Print(saveName);
   
  return 0;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int compareDataMC_VtxDistribution(bool withWeights = true){

  //TeresaPlottingStyle::init();
  gStyle->SetOptStat("");
  
  TFile *file, *file1;  
  char title[100];
  TH1D* hNVtxD[numTrigger];
  TH1D* hNVtxM[numTrigger];
  TCanvas* canv[numTrigger];
  TLatex *info;
  
  for(int i=0; i<numTrigger;i++){

    sprintf(title,"../plots_2012/PF_L1CHS/data/root_files/VtxN%i_PFCHS_data.root",i);
    file = TFile::Open(title);  
    file->GetObject("histo",hNVtxD[i]);
    
    if(withWeights) sprintf(title,"../plots_2012/PF_L1CHS/mc/root_files/VtxN%i_PFCHS_mc.root",i);
    else            sprintf(title,"../plots_2012/PF_L1CHS/mc/root_files/VtxNwoWeights%i_PFCHS_mc.root",i);
    file1 = TFile::Open(title);
    file1->GetObject("histo",hNVtxM[i]);
    
    hNVtxM[i]->Scale(hNVtxD[i]->Integral()/hNVtxM[i]->Integral()); 
    
    if(i == 0)      sprintf(title,"HLT_Photon20_CaloIdVL_IsoL");
    else if(i == 1) sprintf(title,"HLT_Photon30_CaloIdVL_IsoL");
    else if(i == 2) sprintf(title,"HLT_Photon50_CaloIdVL_IsoL");
    else if(i == 3) sprintf(title,"HLT_Photon75_CaloIdVL_IsoL");
    else if(i == 4) sprintf(title,"HLT_Photon90_CaloIdVL_IsoL");
    else if(i == 5) sprintf(title,"HLT_Photon135");
    else if(i == 6) sprintf(title,"HLT_Photon150");

    canv[i] = DrawComparison(hNVtxD[i], hNVtxM[i],title,"","#Vtx", 0);
    canv[i] ->SetLogy();

    info  = new TLatex();
    info -> SetNDC();    
    info -> SetTextSize(0.04);

    if(withWeights) info->DrawLatex(0.65,0.60, "reweighed");
    else       info->DrawLatex(0.65,0.60, "#splitline{not}{reweighed}");

    if(withWeights) sprintf(title,"plots/NVtxComparison%i.pdf",i);
    else            sprintf(title,"plots/NVtxComparisonWoWeights%i.pdf",i);
    canv[i] ->Print(title);
  
  }

  delete file;
  delete file1;
  
  return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int compare_PhotonPtWithWoWeights(){

  TeresaPlottingStyle::init();
  gStyle->SetOptStat("");
  
  TFile *file;  
  TH1D *histo, *histoReweighted;
  TCanvas* canv = new TCanvas("c", "Photon Pt", 700, 700);
  canv->cd();  

  TString filename = "../plots_2012/PF_L1CHS/mc/root_files/hPhoton1PtwoWeights_PFCHS_mc.root";
  file = TFile::Open(filename);
  file->GetObject("histo",histo);
  filename = "../plots_2012/PF_L1CHS/mc/root_files/hPhoton1Pt_PFCHS_mc.root";
  file = TFile::Open(filename);
  file->GetObject("histo",histoReweighted);
      
  histo           -> SetLineColor(9);
  histo           -> SetMarkerColor(9);
  histo           -> SetMarkerStyle(1);
  histo           -> Rebin(10);
  histoReweighted -> SetLineColor(2);
  histoReweighted -> SetMarkerColor(2);
  histoReweighted -> SetMarkerStyle(1);
  histoReweighted -> Rebin(10);

  histo -> SetMinimum(0.000000001);
  histo -> SetMaximum(histo->GetMaximum()*40);
  histo -> SetYTitle("frequency");
  histo -> SetTitle("");

  histo->Draw();
  histoReweighted->Draw("same");

  TLegend* leg1 = new TLegend(0.60, 0.73, 0.90, 0.90);
  leg1->SetFillStyle(0);
  //leg1->SetLineStyle(1);
  //leg1->SetTextFont(42);
  //leg1->SetTextSize(0.05);
  leg1->AddEntry(histo, "not reweighted" , "lep");
  leg1->AddEntry(histoReweighted, "reweighted", "lep");
  leg1->Draw("same");
 

  //canv = DrawComparison(histo, histoReweighted,"photon Pt","","Photon pt", 0);
  canv->SetLogy();
  canv ->Print("plots/PhotonPtComparison_reweighted.pdf");
  
  delete file;
  
  return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Some Functions

TGraphErrors* GetTGraphErrors(TString filename, TString objectName){

  TGraphErrors* graph = 0;
  TFile* file =  TFile::Open(filename);
  if(file != 0) file -> GetObject(objectName,graph);
  delete file;
  
  return graph;
}


TF1* GetTF1(TString filename, TString objectName){

  TF1* object = 0;

  TFile* file =  TFile::Open(filename);
  if(file != 0) file -> GetObject(objectName,object);
  delete file;

  return object;
}

TH1D* GetTH1D(TString filename, TString objectName){

  TH1D* object = 0;

  TFile* file =  TFile::Open(filename);
  if(file != 0) file -> GetObject(objectName,object);
  object -> SetDirectory(0);
  delete file;

  return object;
}

TCanvas* DrawComparison(TH1D* prediction, TH1D* selection, TString Title, TString LumiTitle, TString xTitle, bool withPad2)
{

  double MinX = selection->GetXaxis()->GetXmin();
  double MaxX = selection->GetXaxis()->GetXmax();
  double MaxY = selection->GetMaximum();
  double YRangeMax = 2*MaxY;
  TString titlePrediction;
  TString titleSelection;
  TString RatioTitle;
  

  titlePrediction = "Data";
  titleSelection = "MC";
  RatioTitle = "(Data-MC)/MC";

  //static Int_t c_LightBrown   = TColor::GetColor( "#D9D9CC" );
  static Int_t c_LightGray    = TColor::GetColor( "#DDDDDD" );

  prediction->SetAxisRange(MinX, MaxX, "X");
  prediction->GetYaxis()->SetRangeUser(0.005, YRangeMax);
  prediction->SetMarkerStyle(20);
  prediction->SetMarkerSize(0.9);
  prediction->SetMarkerColor(kBlack);
  prediction->SetXTitle(xTitle);
  prediction->SetYTitle("frequency");
  selection->SetAxisRange(MinX, MaxX, "X");
  selection->GetYaxis()->SetRangeUser(0.05, YRangeMax);
  // selection->SetFillColor(c_LightBrown);
  selection->SetFillColor(c_LightGray);
  selection->SetTitle("");
  selection->SetXTitle(xTitle);
  selection->SetYTitle("frequency");
  selection->SetTitleSize(0.05, "XYZ");
  selection->SetTitleOffset(0.65, "Y");
  selection->SetLabelSize(0.05,"XYZ");
  selection->SetTitleOffset(1.2,"XYZ");
  TCanvas *c = new TCanvas("c", "Comparison and ratio of two histos", 700, 700);
  TPad *pad1;
  if(withPad2) pad1 = new TPad("pad1a", "pad1a", 0, 0.35, 1, 1);
  else         pad1 = new TPad("pad1a", "pad1a", 0, 0, 1, 1);
  pad1->SetLogy();
  if(withPad2) pad1->SetBottomMargin(0);
  else pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);
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
  
  TLegend* leg1 = new TLegend(0.60, 0.73, 0.90, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetLineStyle(1);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.05);
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
  if(withPad2){
    TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 0.35);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.15);
    pad2->Draw();
    pad2->cd();
    TH1D* r = new TH1D(*selection);
    r->SetTitle("");
    r->SetLabelSize(0.1, "XYZ");
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
  }
   return c;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
