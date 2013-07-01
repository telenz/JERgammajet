#include <TROOT.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TFile.h>
#include <TPostScript.h>
#include <TLegend.h>
#include <TMath.h>
#include <TArrayF.h>
#include <TLine.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <vector>
#include <string>

#include "/afs/naf.desy.de/user/t/telenz/comparison/tdrstyle_mod.C"

using namespace std;


void methodComparison2012(){

   setTDRStyle(false);
 
   double eta_bins[5] = {0., 0.5, 1.1, 1.7, 2.3};
   
   TH1F *Res_INPUT = new TH1F("INPUT","", 4, eta_bins);
   TH1F *Res_Gaus  = new TH1F("MCsmeared_MC_ratio_Gaus", "", 4, eta_bins);
   TH1F *Res_RMS99 = new TH1F("MCsmeared_MC_ratio_RMS99", "", 4, eta_bins);
   TH1F *Res_RMS95 = new TH1F("MCsmeared_MC_ratio_RMS95", "", 4, eta_bins);
   TH1F *Res_2011Final = new TH1F("Data_MC_ratio_2011","", 4, eta_bins);


   // 2011
   Res_2011Final->SetBinContent(1, 1.052);
   Res_2011Final->SetBinContent(2, 1.057);
   Res_2011Final->SetBinContent(3, 1.096);
   Res_2011Final->SetBinContent(4, 1.134);
   
   // 2011
   TGraphAsymmErrors *Res_2011 = new TGraphAsymmErrors(Res_2011Final);
   Res_2011->SetPointError(0, 0., 0., 0.063, 0.062);
   Res_2011->SetPointError(1, 0., 0., 0.057, 0.056);
   Res_2011->SetPointError(2, 0., 0., 0.065, 0.064);
   Res_2011->SetPointError(3, 0., 0., 0.094, 0.092);

   //OLD
   /*
   // 2sigma Gauss CHS
   Res_Gaus->SetBinContent(1, 1.134);
   Res_Gaus->SetBinContent(2, 1.096);
   Res_Gaus->SetBinContent(3, 1.158);
   Res_Gaus->SetBinContent(4, 1.287);
   
   Res_Gaus->SetBinError(1, 0.019);
   Res_Gaus->SetBinError(2, 0.019);
   Res_Gaus->SetBinError(3, 0.024);
   Res_Gaus->SetBinError(4, 0.049);


   // RMS 99
   Res_RMS99->SetBinContent(1, 1.174);
   Res_RMS99->SetBinContent(2, 1.187);
   Res_RMS99->SetBinContent(3, 1.202);
   Res_RMS99->SetBinContent(4, 1.419);
   
   Res_RMS99->SetBinError(1, 0.013);
   Res_RMS99->SetBinError(2, 0.014);
   Res_RMS99->SetBinError(3, 0.018);
   Res_RMS99->SetBinError(4, 0.038);


   // RMS 95   
   Res_RMS95->SetBinContent(1, 1.149);
   Res_RMS95->SetBinContent(2, 1.149);
   Res_RMS95->SetBinContent(3, 1.188);
   Res_RMS95->SetBinContent(4, 1.333);
   
   Res_RMS95->SetBinError(1, 0.014);
   Res_RMS95->SetBinError(2, 0.014);
   Res_RMS95->SetBinError(3, 0.017);
   Res_RMS95->SetBinError(4, 0.037);
   */

   // NEW
   // 2sigma Gauss CHS
   Res_Gaus->SetBinContent(1, 1.107);
   Res_Gaus->SetBinContent(2, 1.066);
   Res_Gaus->SetBinContent(3, 1.072);
   Res_Gaus->SetBinContent(4, 1.238);
   
   Res_Gaus->SetBinError(1, 0.021);
   Res_Gaus->SetBinError(2, 0.019);
   Res_Gaus->SetBinError(3, 0.026);
   Res_Gaus->SetBinError(4, 0.043);

   // NEW (with Trigger)
   // 2sigma Gauss CHS
   Res_Gaus->SetBinContent(1, 1.101);
   Res_Gaus->SetBinContent(2, 1.072);
   Res_Gaus->SetBinContent(3, 1.078);
   Res_Gaus->SetBinContent(4, 1.240);
   
   Res_Gaus->SetBinError(1, 0.021);
   Res_Gaus->SetBinError(2, 0.020);
   Res_Gaus->SetBinError(3, 0.028);
   Res_Gaus->SetBinError(4, 0.043);


   // RMS 99
   Res_RMS99->SetBinContent(1, 1.063);
   Res_RMS99->SetBinContent(2, 1.090);
   Res_RMS99->SetBinContent(3, 1.093);
   Res_RMS99->SetBinContent(4, 1.240);
   
   Res_RMS99->SetBinError(1, 0.014);
   Res_RMS99->SetBinError(2, 0.014);
   Res_RMS99->SetBinError(3, 0.019);
   Res_RMS99->SetBinError(4, 0.041);

   // RMS 99 (with Trigger in MC)
   Res_RMS99->SetBinContent(1, 1.063);
   Res_RMS99->SetBinContent(2, 1.091);
   Res_RMS99->SetBinContent(3, 1.094);
   Res_RMS99->SetBinContent(4, 1.241);
   
   Res_RMS99->SetBinError(1, 0.014);
   Res_RMS99->SetBinError(2, 0.014);
   Res_RMS99->SetBinError(3, 0.014);
   Res_RMS99->SetBinError(4, 0.041);


   // RMS 95   
   Res_RMS95->SetBinContent(1, 1.068);
   Res_RMS95->SetBinContent(2, 1.082);
   Res_RMS95->SetBinContent(3, 1.114);
   Res_RMS95->SetBinContent(4, 1.238);
   
   Res_RMS95->SetBinError(1, 0.012);
   Res_RMS95->SetBinError(2, 0.015);
   Res_RMS95->SetBinError(3, 0.018);
   Res_RMS95->SetBinError(4, 0.040);

   // RMS 95   (with Trigger in MC)
   Res_RMS95->SetBinContent(1, 1.063);
   Res_RMS95->SetBinContent(2, 1.083);
   Res_RMS95->SetBinContent(3, 1.117);
   Res_RMS95->SetBinContent(4, 1.238);
   
   Res_RMS95->SetBinError(1, 0.012);
   Res_RMS95->SetBinError(2, 0.015);
   Res_RMS95->SetBinError(3, 0.018);
   Res_RMS95->SetBinError(4, 0.040);

   // RMS 95   With Pixel Seed 
   //Res_RMS95->SetBinContent(1, 1.068);
   //Res_RMS95->SetBinContent(2, 1.082);
   //Res_RMS95->SetBinContent(3, 1.114);
   //Res_RMS95->SetBinContent(4, 1.238);
   
   //Res_RMS95->SetBinError(1, 0.012);
   //Res_RMS95->SetBinError(2, 0.015);
   //Res_RMS95->SetBinError(3, 0.018);
   //Res_RMS95->SetBinError(4, 0.040);

   //Without PixelSeed RMS 95
   /*
   Res_Gaus->SetBinContent(1, 1.145);
   Res_Gaus->SetBinContent(2, 1.153);
   Res_Gaus->SetBinContent(3, 1.172);
   Res_Gaus->SetBinContent(4, 1.336);
   
   Res_Gaus->SetBinError(1, 0.013);
   Res_Gaus->SetBinError(2, 0.011);
   Res_Gaus->SetBinError(3, 0.018);
   Res_Gaus->SetBinError(4, 0.037);
   */

   ///////////////////////////////////////////////////////////////
  
   TCanvas *c = new TCanvas();
   Res_Gaus->SetXTitle("|#eta|");
   Res_Gaus->GetXaxis()->SetRangeUser(0., 5.);
   Res_Gaus->SetYTitle("Data/MC ratio (const fit)");
   Res_Gaus->GetYaxis()->SetRangeUser(0.8, 1.5);
   Res_Gaus->SetMarkerStyle(21);
   //Res_Gaus->SetMarkerStyle(25);
   Res_Gaus->SetMarkerSize(1.2);
   //Res_Gaus->SetLineColor(kBlack);
   //Res_Gaus->SetMarkerColor(kBlack);
   Res_Gaus->SetLineColor(kRed);
   Res_Gaus->SetMarkerColor(kRed);
   Res_Gaus->Draw("e1p");
   Res_2011->SetMarkerStyle(20);
   Res_2011->SetMarkerSize(1.4);
   Res_2011->SetFillColor(kGray);
   Res_2011->SetFillStyle(3001);
   Res_2011->SetLineColor(kGray);
   Res_2011->DrawClone("e3psame");


   Res_RMS99->SetMarkerStyle(23);
   Res_RMS99->SetMarkerSize(1.2);
   Res_RMS99->SetLineColor(kBlue-4);
   Res_RMS99->SetMarkerColor(kBlue-4);
   Res_RMS99->Draw("e1psame");

   Res_RMS95->SetMarkerStyle(24);
   Res_RMS95->SetMarkerSize(1.2);
   Res_RMS95->SetLineColor(3);
   Res_RMS95->SetMarkerColor(3);
   Res_RMS95->Draw("e1psame");

   Res_Gaus->Draw("e1psame");

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

   leg->AddEntry(Res_2011Final,"2011 (total error)", "pfl");
   leg->AddEntry(Res_Gaus,"Gaus Fit 2 #sigma", "pfl");
   leg->AddEntry(Res_RMS99,"RMS 99%", "pfl");
   leg->AddEntry(Res_RMS95,"RMS 95%", "pfl");

   //leg->AddEntry(Res_Gaus,"RMS 95% without Pixel Seed Veto", "pfl");
   //leg->AddEntry(Res_RMS95,"RMS 95% with Pixel Seed Veto", "pfl");

   leg->Draw("same");
   c->Print("plots/resultsComparison.pdf");   
   //c->Print("MCClosure.pdf");   
}
