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


void PlotResolutions_2011(){
   setTDRStyle();

   tdrStyle->SetNdivisions(505, "X");

 
   double eta_bins[5] = {0., 0.5, 1.1, 1.7, 2.3};
   
   TH1F *Res_2011Final = new TH1F("Data_MC_ratio_2011","", 4, eta_bins);
   TH1F *Res_2012AB = new TH1F("Data_MC_ratio_2012AB", "", 4, eta_bins);
   TH1F *Res_2012ABC = new TH1F("Data_MC_ratio_2012ABC", "", 4, eta_bins); // prompt reco datasets
   TH1F *Res_2012ABC_CHS = new TH1F("Data_MC_ratio_2012ABC_CHS", "", 4, eta_bins); // prompt reco datasets

   TH1F *Res_2012ABC_lowPU = new TH1F("Data_MC_ratio_2012ABC_lowPU", "", 4, eta_bins);
   TH1F *Res_2012ABC_highPU = new TH1F("Data_MC_ratio_2012ABC_highPU", "", 4, eta_bins);

   TH1F *Res_2012ABC_rereco_2012SQLV7 = new TH1F("Data_MC_ratio_2012ABC_rereco_2012SQLV7", "", 4, eta_bins); // A-C1 rereco, C2 prompt reco (HCP)
   TH1F *Res_2012ABC_rereco_2012FallV4 = new TH1F("Data_MC_ratio_2012ABC_rereco_2012FallV4", "", 4, eta_bins);
   TH1F *Res_2012ABC_rereco_2012SQLV7_chs = new TH1F("Data_MC_ratio_2012ABC_rereco_2012SQLV7_chs", "", 4, eta_bins); 
   TH1F *Res_2012ABC_rereco_2012FallV4_chs = new TH1F("Data_MC_ratio_2012ABC_rereco_2012FallV4_chs", "", 4, eta_bins);

   Res_2011Final->SetBinContent(1, 1.052);
   Res_2011Final->SetBinContent(2, 1.057);
   Res_2011Final->SetBinContent(3, 1.096);
   Res_2011Final->SetBinContent(4, 1.134);
   

   TGraphAsymmErrors *Res_2011 = new TGraphAsymmErrors(Res_2011Final);
   Res_2011->SetPointError(0, 0., 0., 0.063, 0.062);
   Res_2011->SetPointError(1, 0., 0., 0.057, 0.056);
   Res_2011->SetPointError(2, 0., 0., 0.065, 0.064);
   Res_2011->SetPointError(3, 0., 0., 0.094, 0.092);
   

   Res_2012AB->SetBinContent(1, 1.102);
   Res_2012AB->SetBinContent(2, 1.070);
   Res_2012AB->SetBinContent(3, 1.098);
   Res_2012AB->SetBinContent(4, 1.141);
   
   Res_2012AB->SetBinError(1, 0.011);
   Res_2012AB->SetBinError(2, 0.012);
   Res_2012AB->SetBinError(3, 0.020);
   Res_2012AB->SetBinError(4, 0.042);
   
   //Old data
   Res_2012ABC->SetBinContent(1, 1.072);
   Res_2012ABC->SetBinContent(2, 1.080);
   Res_2012ABC->SetBinContent(3, 1.038);
   Res_2012ABC->SetBinContent(4, 1.229);

   //Old data Errors
   Res_2012ABC->SetBinError(1, 0.025);
   Res_2012ABC->SetBinError(2, 0.025);
   Res_2012ABC->SetBinError(3, 0.030);
   Res_2012ABC->SetBinError(4, 0.068);

   //New data
   Res_2012ABC->SetBinContent(1, 1.056);
   Res_2012ABC->SetBinContent(2, 1.092);
   Res_2012ABC->SetBinContent(3, 1.147);
   Res_2012ABC->SetBinContent(4, 1.165);
   
   //New data Errors
   Res_2012ABC->SetBinError(1, 0.024);
   Res_2012ABC->SetBinError(2, 0.024);
   Res_2012ABC->SetBinError(3, 0.030);
   Res_2012ABC->SetBinError(4, 0.152);

   
   // Old data
   Res_2012ABC_CHS->SetBinContent(1, 1.031);
   Res_2012ABC_CHS->SetBinContent(2, 1.060);
   Res_2012ABC_CHS->SetBinContent(3, 1.087);
   Res_2012ABC_CHS->SetBinContent(4, 1.216);
   
   // Old data errors
   Res_2012ABC_CHS->SetBinError(1, 0.021);
   Res_2012ABC_CHS->SetBinError(2, 0.028);
   Res_2012ABC_CHS->SetBinError(3, 0.023);
   Res_2012ABC_CHS->SetBinError(4, 0.050);
   
   // New data
   Res_2012ABC_CHS->SetBinContent(1, 1.144);
   Res_2012ABC_CHS->SetBinContent(2, 1.160);
   Res_2012ABC_CHS->SetBinContent(3, 1.197);
   Res_2012ABC_CHS->SetBinContent(4, 1.340);
   
   // New data errors
   Res_2012ABC_CHS->SetBinError(1, 0.012);
   Res_2012ABC_CHS->SetBinError(2, 0.011);
   Res_2012ABC_CHS->SetBinError(3, 0.015);
   Res_2012ABC_CHS->SetBinError(4, 0.034);

   
   // New data 2011
   //Res_2012ABC_CHS->SetBinContent(1, 1.056);
   //Res_2012ABC_CHS->SetBinContent(2, 1.092);
   //Res_2012ABC_CHS->SetBinContent(3, 1.147);
   //Res_2012ABC_CHS->SetBinContent(4, 1.165);
   
   // New data errors 2011
   //Res_2012ABC_CHS->SetBinError(1, 0.024);
   //Res_2012ABC_CHS->SetBinError(2, 0.024);
   //Res_2012ABC_CHS->SetBinError(3, 0.030);
   //Res_2012ABC_CHS->SetBinError(4, 0.152);
   

   Res_2012ABC_lowPU->SetBinContent(1, 1.105);
   Res_2012ABC_lowPU->SetBinContent(2, 1.102);
   Res_2012ABC_lowPU->SetBinContent(3, 1.165);
   Res_2012ABC_lowPU->SetBinContent(4, 1.245);
   
   Res_2012ABC_lowPU->SetBinError(1, 0.009);
   Res_2012ABC_lowPU->SetBinError(2, 0.010);
   Res_2012ABC_lowPU->SetBinError(3, 0.016);
   Res_2012ABC_lowPU->SetBinError(4, 0.035);
   

   Res_2012ABC_highPU->SetBinContent(1, 1.086);
   Res_2012ABC_highPU->SetBinContent(2, 1.075);
   Res_2012ABC_highPU->SetBinContent(3, 1.134);
   Res_2012ABC_highPU->SetBinContent(4, 1.224);
   
   Res_2012ABC_highPU->SetBinError(1, 0.011);
   Res_2012ABC_highPU->SetBinError(2, 0.013);
   Res_2012ABC_highPU->SetBinError(3, 0.021);
   Res_2012ABC_highPU->SetBinError(4, 0.048);
   
   Res_2012ABC_rereco_2012SQLV7->SetBinContent(1, 1.087);
   Res_2012ABC_rereco_2012SQLV7->SetBinContent(2, 1.111);
   Res_2012ABC_rereco_2012SQLV7->SetBinContent(3, 1.170);
   Res_2012ABC_rereco_2012SQLV7->SetBinContent(4, 1.307);
   
   Res_2012ABC_rereco_2012SQLV7->SetBinError(1, 0.007);
   Res_2012ABC_rereco_2012SQLV7->SetBinError(2, 0.007);
   Res_2012ABC_rereco_2012SQLV7->SetBinError(3, 0.013);
   Res_2012ABC_rereco_2012SQLV7->SetBinError(4, 0.028);
   

   Res_2012ABC_rereco_2012FallV4->SetBinContent(1, 1.093);
   Res_2012ABC_rereco_2012FallV4->SetBinContent(2, 1.104);
   Res_2012ABC_rereco_2012FallV4->SetBinContent(3, 1.154);
   Res_2012ABC_rereco_2012FallV4->SetBinContent(4, 1.285);
   
   Res_2012ABC_rereco_2012FallV4->SetBinError(1, 0.007);
   Res_2012ABC_rereco_2012FallV4->SetBinError(2, 0.007);
   Res_2012ABC_rereco_2012FallV4->SetBinError(3, 0.012);
   Res_2012ABC_rereco_2012FallV4->SetBinError(4, 0.028);
   

   Res_2012ABC_rereco_2012SQLV7_chs->SetBinContent(1, 1.092);
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinContent(2, 1.109);
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinContent(3, 1.162);
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinContent(4, 1.259);
   
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinError(1, 0.007);
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinError(2, 0.007);
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinError(3, 0.012);
   Res_2012ABC_rereco_2012SQLV7_chs->SetBinError(4, 0.026);
   

   Res_2012ABC_rereco_2012FallV4_chs->SetBinContent(1, 1.094);
   Res_2012ABC_rereco_2012FallV4_chs->SetBinContent(2, 1.109);
   Res_2012ABC_rereco_2012FallV4_chs->SetBinContent(3, 1.157);
   Res_2012ABC_rereco_2012FallV4_chs->SetBinContent(4, 1.300);
   
   Res_2012ABC_rereco_2012FallV4_chs->SetBinError(1, 0.007);
   Res_2012ABC_rereco_2012FallV4_chs->SetBinError(2, 0.007);
   Res_2012ABC_rereco_2012FallV4_chs->SetBinError(3, 0.011);
   Res_2012ABC_rereco_2012FallV4_chs->SetBinError(4, 0.026);
   

   ///////////////////////////////////////////////////////////////
  
   TCanvas *c = new TCanvas();
   Res_2012ABC_CHS->SetXTitle("|#eta|");
   Res_2012ABC_CHS->GetXaxis()->SetRangeUser(0., 5.);
   Res_2012ABC_CHS->SetYTitle("Data/MC ratio (const fit)");
   Res_2012ABC_CHS->GetYaxis()->SetRangeUser(0.8, 1.5);
   Res_2012ABC_CHS->SetMarkerStyle(21);
   Res_2012ABC_CHS->SetMarkerSize(1.4);
   Res_2012ABC_CHS->SetLineColor(kPink-3);
   Res_2012ABC_CHS->SetMarkerColor(kPink-3);
   Res_2012ABC_CHS->Draw("e1p");
   Res_2011->SetMarkerStyle(20);
   Res_2011->SetMarkerSize(1.4);
   Res_2011->SetFillColor(kGray);
   Res_2011->SetFillStyle(3001);
   Res_2011->SetLineColor(kGray);
   Res_2011->DrawClone("e3psame");
   Res_2012ABC_CHS->Draw("e1psame");
   Res_2012ABC->SetMarkerStyle(23);
   Res_2012ABC->SetMarkerSize(1.4);
   //  Res_2012ABC->SetMarkerColor(kBlue-3);
   Res_2012ABC->SetLineColor(kBlue-4);
   Res_2012ABC->SetMarkerColor(kBlue-4);
   //Res_2012ABC->Draw("e1psame");

   Res_2011->SetPointError(0, 0., 0., 0., 0.);
   Res_2011->SetPointError(1, 0., 0., 0., 0.);
   Res_2011->SetPointError(2, 0., 0., 0., 0.);
   Res_2011->SetPointError(3, 0., 0., 0., 0.);
   Res_2011->SetPointError(4, 0., 0., 0., 0.);

   Res_2011->Draw("psame");

   cmsPrel();

   TLegend *leg = new TLegend(0.54, 0.17, 0.74, 0.4);
   leg->SetBorderSize(0);
   // leg->SetBorderMode(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.035);

   leg->AddEntry(Res_2011,"2011 (total error)", "pfl");
   //leg->AddEntry(Res_2012ABC_CHS,"#splitline{2012 ABC CHSJets}{(stat. error)}", "P");
   leg->AddEntry(Res_2012ABC_CHS,"2011 #gamma + Jet RMS 95%", "P");
   //leg->AddEntry(Res_2012ABC_CHS,"2011 #gamma + Jet", "P");
   //leg->AddEntry(Res_2012ABC,"2012 ABC (stat. error)", "P");
   //leg->AddEntry(Res_2012ABC,"2011 #gamma + Jet 2#sigma Gaus", "P");

   leg->Draw("same");

   
   c->Print("CompareResolutions_total.pdf");


   
}
