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
#include "../CODE/myFunctions.h"
#include "../CODE/myFunctions.C"
#include "plotStyle.h"
#include "/afs/naf.desy.de/user/t/telenz/comparison/tdrstyle_mod.C"


int plotDataMCComparison(int eta){


  TeresaPlottingStyle::init();

  const TString method = "RMS99";

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
    Ratio   -> Fit("name","QR0");
   
   
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
    band->SetLineColor(kRed-9);
    band->SetFillColor(kRed-9);
    band->SetMarkerColor(kRed-9);
    band->SetFillStyle(1001);
    band->SetMarkerSize(0.00001);
    band ->SetMinimum(0.5);
    band ->SetMaximum(0.5);
    band ->Draw("e3same");
    Ratio  -> Draw("Pesame"); 
  
    TLatex*  info   = new TLatex();
    info -> SetTextFont(132);
    info -> SetNDC();
    info -> SetTextSize(0.05); 
    info -> DrawLatex(0.20,0.78,Form("#chi^{2}/ndof = %4.2f / %i",f1 -> GetChisquare(),f1 -> GetNDF()));
    info -> DrawLatex(0.60,0.78,Form("f = %4.3f #pm %4.3f", f1 -> GetParameter(0), f1->GetParError(0)));
    info -> DrawLatex(0.60,0.85,Form("#bf{%4.1f < |#eta^{Jet}| < %4.1f}",etaBins[eta],etaBins[eta+1]));
    double m_lumi = 19.7;

    TString infotext = TString::Format("CMS Preliminary, %3.1f fb^{-1}", m_lumi);
    TLatex *text1 = new TLatex(3.5, 24, infotext);
    text1->SetNDC();
    text1->SetX(0.22);
    text1->SetTextFont(42);
    infotext = TString::Format("#sqrt{s} = 8 TeV");
    TLatex *text2 = new TLatex(3.5, 24, infotext);
    text2->SetNDC();
    text2->SetX(0.22);
    text2->SetTextFont(42);

    text1->SetTextSize(0.040);
    text1->SetTextAlign(11);
    text1->SetY(0.96);
    text1->SetX(0.15);

    text2->SetTextSize(0.040);
    text2->SetTextAlign(31);
    text2->SetY(0.96);
    text2->SetX(0.95);

    text1->Draw("same");
    text2->Draw("same");

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
