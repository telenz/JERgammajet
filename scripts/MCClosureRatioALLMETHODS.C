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


void MCClosureRatioALLMETHODS(){

   setTDRStyle(true);
   //tdrStyle->SetNdivisions(505, "X");
 
   double eta_bins[5] = {0., 0.5, 1.1, 1.7, 2.3};
   
   TH1F *Res_INPUT = new TH1F("INPUT","", 4, eta_bins);
   TH1F *Res_Gaus  = new TH1F("MCsmeared_MC_ratio_Gaus", "", 4, eta_bins);
   TH1F *Res_RMS99 = new TH1F("MCsmeared_MC_ratio_RMS99", "", 4, eta_bins);
   TH1F *Res_RMS95 = new TH1F("MCsmeared_MC_ratio_RMS95", "", 4, eta_bins);

   // Input Parameters
   Res_INPUT->SetBinContent(1, 1.133);
   Res_INPUT->SetBinContent(2, 1.083);
   Res_INPUT->SetBinContent(3, 1.145);
   Res_INPUT->SetBinContent(4, 1.288);

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

   // RMS 99
   Res_RMS99->SetBinContent(1, 1.048);
   Res_RMS99->SetBinContent(2, 1.068);
   Res_RMS99->SetBinContent(3, 1.086);
   Res_RMS99->SetBinContent(4, 1.116);
   
   Res_RMS99->SetBinError(1, 0.004);
   Res_RMS99->SetBinError(2, 0.004);
   Res_RMS99->SetBinError(3, 0.005);
   Res_RMS99->SetBinError(4, 0.013);

   // RMS 95   
   Res_RMS95->SetBinContent(1, 1.048);
   Res_RMS95->SetBinContent(2, 1.067);
   Res_RMS95->SetBinContent(3, 1.085);
   Res_RMS95->SetBinContent(4, 1.106);
   
   Res_RMS95->SetBinError(1, 0.004);
   Res_RMS95->SetBinError(2, 0.004);
   Res_RMS95->SetBinError(3, 0.006);
   Res_RMS95->SetBinError(4, 0.013);
   

   ///////////////////////////////////////////////////////////////
  
   TCanvas *c = new TCanvas();
   Res_INPUT->SetXTitle("|#eta|");
   Res_INPUT->GetXaxis()->SetRangeUser(0., 5.);
   Res_INPUT->SetYTitle("MC_{smeared}/MC ratio (const fit)");
   Res_INPUT->GetYaxis()->SetRangeUser(1.04, 1.14);
   Res_INPUT->SetMarkerStyle(20);
   Res_INPUT->SetMarkerSize(1.2);
   Res_INPUT->SetLineColor(kBlack);
   Res_INPUT->SetMarkerColor(kBlack);
   Res_INPUT->Draw("e1p");

   Res_Gaus->SetMarkerStyle(21);
   Res_Gaus->SetMarkerSize(1.2);
   Res_Gaus->SetLineColor(kPink-3);
   Res_Gaus->SetMarkerColor(kPink-3);
   Res_Gaus->DrawClone("e1psame");

   Res_INPUT->Draw("e1psame");

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

   cmsPrel();

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
   c->Print("plots/MCClosure.pdf");   
   //c->Print("MCClosure.pdf");   
}
