// Evaluate MC (Out-of-Cone showering) Uncertainty
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

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
#include "../../../../CODE/myDeclarations.h"
#include "../../../../scripts/plotStyle.h"
#include "../../../../CODE/myFunctions.h"
#include "../../../../CODE/myFunctions.C"


int evalMCUncertainty(){

  cout<<endl<<endl<<endl<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for MC uncertainty is executed! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl;
  gErrorIgnoreLevel = 1001;

  TeresaPlottingStyle::init();

  const TString type   = "PFCHS";
  const TString method = "RMS99";

  TString tot_filename, AuxString, fitName;

  TString etaString = Form("Jet Energy Resolution ak7PFCHS - ak5PFCHS");

  double finalErrors[nEtaBins] = {0};
  double finalErrorsE[nEtaBins] = {0};

  TCanvas *c   = 0 ;
  TCanvas *c11 = 0;

  TString pathNameAk7 = "root_files_MCUncertainty/ak7PFCHS_";
  TString pathNameAk5 = "../root_files_FINAL_";

  //pathNameAk7 = "root_files_MCUncertainty/ak7PFCHS_";
  //pathNameAk5 = "root_files_MCUncertainty/ak5PFCHS_";

  cout<<"root files from following folders:"<<endl<<pathNameAk5<<endl<<pathNameAk7<<endl<<endl<<endl;

  for(int eta=0; eta<nEtaBins;eta++){

    TString rootFileMCak5PFCHS   = pathNameAk5 + (TString) "mc/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + type + (TString) "_mc_" + method + (TString)".root";
    TString rootFileDataak5PFCHS = pathNameAk5 + (TString) "data/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + type + (TString) "_data_" + method + (TString)".root"; 

    TMultiGraph* mg = new TMultiGraph();
    mg -> SetTitle(etaString);
    

    c = new TCanvas("c",etaString,200,10,800,800);
    c -> cd();
    c -> SetBottomMargin(0.14);
    c -> SetLeftMargin(0.14);

    double maximum = 0;
   
    TGraphErrors* graph[4]; 
  
    graph[0] = readTGraphErrors(rootFileMCak5PFCHS,"Graph","Graph");    
    graph[1] = readTGraphErrors(rootFileDataak5PFCHS,"Graph","Graph");   
    graph[2] = readTGraphErrors(rootFileMCak7PFCHS,"Graph","Graph");    
    graph[3] = readTGraphErrors(rootFileDataak7PFCHS,"Graph","Graph");    


    for(int i=0; i<4 ; i++){
      if(graph[i]->GetYaxis()->GetXmax() > maximum) maximum = graph[i]->GetYaxis()->GetXmax();
    }
  
    graph[0] -> SetMarkerColor(8);
    graph[1] -> SetMarkerColor(9);
    graph[2] -> SetMarkerColor(1);
    graph[3] -> SetMarkerColor(2);
 
    graph[0] -> SetLineColor(8);
    graph[1] -> SetLineColor(9);
    graph[2] -> SetLineColor(1);
    graph[3] -> SetLineColor(2);

    graph[0] -> GetFunction("fResolution") -> SetLineColor(8);
    graph[1] -> GetFunction("fResolution") -> SetLineColor(9);
    graph[2] -> GetFunction("fResolution") -> SetLineColor(1);
    graph[3] -> GetFunction("fResolution") -> SetLineColor(2);
  
 
    TLegend *legend  = new TLegend(0.5,0.7,0.9,0.9);
    legend -> SetFillColor(0);
    legend -> SetTextSize(0.033);
    legend -> AddEntry(graph[0],"MC ak5 PFCHS","l");
    legend -> AddEntry(graph[2],"MC ak7 PFCHS","l");
    legend -> AddEntry(graph[1],"Data ak5 PFCHS","l");
    legend -> AddEntry(graph[3],"Data ak7 PFCHS","l");
  
    mg->Add(graph[0]);  
    mg->Add(graph[1]);  
    mg->Add(graph[2]); 
    mg->Add(graph[3]);   
  
    mg->Draw("AP");
    
    mg -> GetYaxis() -> SetTitle("JER");
    mg -> SetMinimum(0.02);
    mg -> SetMaximum(0.12);   
    mg -> GetXaxis() -> SetLimits(0,600);
    mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]"); 
    legend->Draw("same");
  
    TLatex*  info   = new TLatex();
    info-> SetNDC();
    AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta],etaBins[eta+1]);
    info->DrawLatex(0.6,0.25,AuxString);

    tot_filename = (TString) "plotsMC/Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_MCUncertainty_" + method + (TString)".pdf";
    c -> Print(tot_filename);
    delete c;
  
    // -------------------  Ratio Plot for Resolution (Calculation for ak5 and ak7 PFCHS)  -------------------
    double xRatioak5[nPtBins]  = {0};
    double yRatioak5[nPtBins]  = {0};
    double xERatioak5[nPtBins] = {0};
    double yERatioak5[nPtBins] = {0};

    double xRatioak7[nPtBins]  = {0};
    double yRatioak7[nPtBins]  = {0};
    double xERatioak7[nPtBins] = {0};
    double yERatioak7[nPtBins] = {0};


    double* ak5dataX  = graph[1] -> GetX();
    double* ak5dataY  = graph[1] -> GetY();
    double* ak5dataEX = graph[1] -> GetEX();
    double* ak5dataEY = graph[1] -> GetEY();

    double* ak5mcX  = graph[0] -> GetX();
    double* ak5mcY  = graph[0] -> GetY();
    double* ak5mcEX = graph[0] -> GetEX();
    double* ak5mcEY = graph[0] -> GetEY();

    double* ak7dataX  = graph[3] -> GetX();
    double* ak7dataY  = graph[3] -> GetY();
    double* ak7dataEX = graph[3] -> GetEX();
    double* ak7dataEY = graph[3] -> GetEY();

    double* ak7mcX  = graph[2] -> GetX();
    double* ak7mcY  = graph[2] -> GetY();
    double* ak7mcEX = graph[2] -> GetEX();
    double* ak7mcEY = graph[2] -> GetEY();

  
    int idxMC = 0;
    for(int i=0; i<graph[1]->GetN(); i++){
  
      if(abs(ak5dataX[i]/ak5mcX[idxMC] -1.) > 0.05){
	idxMC += 1;
	i     -= 1;
	continue;
      }
    
      xRatioak5[i]  = 1./2.*(ak5dataX[i] + ak5mcX[idxMC]);
      yRatioak5[i]  = ak5dataY[i]/ak5mcY[idxMC];
      xERatioak5[i] =  1./2.*TMath::Sqrt(TMath::Power(ak5dataEX[i],2)+TMath::Power(ak5mcEX[idxMC],2));
      yERatioak5[i] = sqrt(TMath::Power((1./ak5mcY[idxMC]),2)*TMath::Power(ak5dataEY[i],2)+TMath::Power((ak5dataY[i]/(TMath::Power(ak5mcY[idxMC],2))),2)*TMath::Power(ak5mcEY[idxMC],2));
   
      idxMC += 1;
    }


    idxMC = 0;
    for(int i=0; i<graph[3]->GetN(); i++){
  
      if(abs(ak7dataX[i]/ak7mcX[idxMC] -1.) > 0.05){
	idxMC += 1;
	i     -= 1;
	continue;
      }

      xRatioak7[i]  = 1./2.*(ak7dataX[i] + ak7mcX[idxMC]);
      yRatioak7[i]  = ak7dataY[i]/ak7mcY[idxMC];
      xERatioak7[i] =  1./2.*TMath::Sqrt(TMath::Power(ak7dataEX[i],2)+TMath::Power(ak7mcEX[idxMC],2));
      yERatioak7[i] = sqrt(TMath::Power((1./ak7mcY[idxMC]),2)*TMath::Power(ak7dataEY[i],2)+TMath::Power((ak7dataY[i]/(TMath::Power(ak7mcY[idxMC],2))),2)*TMath::Power(ak7mcEY[idxMC],2));
   
      idxMC += 1;
    }


    TGraphErrors *Ratioak5 = new TGraphErrors(graph[1]->GetN(),xRatioak5,yRatioak5,xERatioak5,yERatioak5);
    TGraphErrors *Ratioak7 = new TGraphErrors(graph[3]->GetN(),xRatioak7,yRatioak7,xERatioak7,yERatioak7);

    TString ratioName = "Data/MC Ratio";
    Ratioak5 -> SetTitle(ratioName); 
    ratioName = "Data/MC Ratio - ak7 PFCHS";
    Ratioak7 -> SetTitle(ratioName); 
    Ratioak5 -> GetXaxis()->SetLimits(0,600);
    Ratioak7 -> GetXaxis()->SetLimits(0,600); 
    Ratioak5 -> GetXaxis() -> SetTitle("Photon pT");
    Ratioak5 -> GetYaxis() -> SetTitle("Ratio of JER (Data/MC)");
    Ratioak7 -> GetXaxis() -> SetTitle("Photon pT");
    Ratioak7 -> GetYaxis() -> SetTitle("Ratio of JER (Data/MC)");
    Ratioak5 -> SetMinimum(0.8);
    Ratioak5 -> SetMaximum(1.4);
    Ratioak7 -> SetMinimum(0.8);
    Ratioak7 -> SetMaximum(1.4);

    TF1* f1ak5 = new TF1("nameak5","pol0",0,600);   
    Ratioak5 -> Fit("nameak5","Q");
    //cout<<"ChiSquare/NDF (ak5) = "<<f1ak5 -> GetChisquare()<<" / "<<f1ak5->GetNDF()<<endl;
    TF1* f1ak7 = new TF1("nameak7","pol0",0,600);   
    Ratioak7 -> Fit("nameak7","Q");
    //cout<<"ChiSquare/NDF (ak7) = "<<f1ak7 -> GetChisquare()<<" / "<<f1ak7->GetNDF()<<endl<<endl;


    Ratioak5 -> SetMarkerColor(1);
    Ratioak5 -> SetMarkerStyle(20);
    Ratioak5 -> SetLineColor(1);
    Ratioak5 -> GetFunction("nameak5")->SetLineColor(1);

    Ratioak7 -> SetMarkerColor(9);
    Ratioak7 -> SetMarkerStyle(24);
    Ratioak7 -> SetLineColor(9);
    Ratioak7 -> GetFunction("nameak7")->SetLineColor(9);

    legend  = 0;
    legend = new TLegend(0.3,0.8,0.9,0.9);
    legend -> SetFillColor(0);
    char legname[100]; 
    TString ak5Entry = (TString) "ak5 PFCHS: " + Form("%4.3f #pm %4.3f", f1ak5 -> GetParameter(0), f1ak5->GetParError(0));
    TString ak7Entry = (TString) "ak7 PFCHS: " + Form("%4.3f #pm %4.3f", f1ak7 -> GetParameter(0), f1ak7->GetParError(0));
    legend->AddEntry(Ratioak5,ak5Entry,"pl");
    legend->AddEntry(Ratioak7,ak7Entry,"pl");
    c11 = new TCanvas("c11",ratioName,800,10,700,700);
    c11 -> cd();
  
    // Draw info boxes
    Ratioak5  -> Draw("AP"); 
    Ratioak7  -> Draw("sameP");
    legend    -> Draw("same");

    TLatex*  info1  = new TLatex();  
    info1 -> SetNDC();
    sprintf(legname,"#splitline{#chi^{2} = %4.2f}{dof = %i}",f1ak5 -> GetChisquare(),f1ak5 -> GetNDF());

    AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBins[eta],etaBins[eta+1]);
    info1->DrawLatex(0.6,0.7,AuxString);

    double ak7      = f1ak7 -> GetParameter(0);
    double ak7Error = f1ak7 -> GetParError(0);
    double ak5      = f1ak5 -> GetParameter(0);
    double ak5Error = f1ak5 -> GetParError(0);

    AuxString = Form("rel. Difference = %4.3f", abs(ak7/ak5-1.));
    info1->DrawLatex(0.4,0.3,AuxString);

    TString pdfFile = (TString) "plotsMC/Ratio_Resolution_for_" + (long) (eta+1) + (TString) "_eta_bin_" + type + (TString) "_" + method + (TString) "_ak5_ak7_comparison.pdf";
    c11 -> SaveAs(pdfFile);
    delete c11;  

    finalErrors[eta]  = abs(ak7/ak5-1.);
    finalErrorsE[eta] = TMath::Sqrt(TMath::Power(1./ak5,2)*TMath::Power(ak7Error,2) + TMath::Power(ak7/TMath::Power(ak5,2),2)*TMath::Power(ak5Error,2) - 2.0 * ak7/TMath::Power(ak5,3)*ak7Error*ak5Error);

    cout<<"Error in "<<eta+1<<". Bin :  "<<fixed<<setprecision(3)<<finalErrors[eta]<<" +/- "<<finalErrorsE[eta]<<endl<<endl;

  }
  
  // Save relative uncertainties for every eta bin in another root-file
  double eta[nEtaBins] = {1.};
  double etaError[nEtaBins] = {0.};
  for(int i =0; i<nEtaBins-1; i++) eta[i+1] = eta[i]+1.;

  
  TGraphErrors* finalErrorsMC = new TGraphErrors(nEtaBins,eta,finalErrors,etaError,finalErrorsE);
  finalErrorsMC -> SetMarkerStyle(20);
  finalErrorsMC -> SetTitle("Final relative Erros on ratio (MC)");
  finalErrorsMC -> GetXaxis() -> SetTitle("#eta^{Jet}");
  finalErrorsMC -> GetYaxis() -> SetTitle("|JER_{ak7} /JER_{ak5} - 1|");
  finalErrorsMC -> GetXaxis() -> SetTitleSize(0.05);
  finalErrorsMC -> GetYaxis() -> SetTitleSize(0.05);  


  tot_filename = "plotsMC/FinalErrosMC.root"; 
  tot_filename = (TString) "plotsMC/FinalErrorsMC_" + type + (TString) "_" + method + (TString) ".root"; 
  TFile *f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(finalErrorsMC,"graph");
  f->Close();
  delete f;

   
  return 0;
}
