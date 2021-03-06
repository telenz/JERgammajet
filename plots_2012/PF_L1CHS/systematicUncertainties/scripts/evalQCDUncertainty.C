// ==== Evaluation of QCD Uncertainty =====
// ========================================

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
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "../../../../CODE/myDeclarations.h"
#include "../../../../scripts/plotStyle.h"
#include "../../../../CODE/myFunctions.h"
#include "../../../../CODE/myFunctions.C"


int evalQCDUncertainty(int eta){

  cout<<endl<<endl<<endl<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for QCD uncertainty is executed! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl<<endl;
  gErrorIgnoreLevel = 1001;
  
  TeresaPlottingStyle::init();
  
  const TString method = "RMS99";
  const TString type   = "PFCHS";
  
  const double etaBinsQCD[3] = {0.0,1.3,2.7};

  TString pathName    = "root_files_QCDUncertainty/mcWoQCD/";
  TString pathNameQCD = "root_files_QCDUncertainty/dataWithQCD/";

  cout<<"root files from following folders:"<<endl<<pathName<<endl<<pathNameQCD<<endl<<endl<<endl;


  TString etaString,tot_filename, AuxString;
  TMultiGraph* mg = new TMultiGraph();

  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Calculate Correlation with the help of the control plot hPhoton1Pt_PFCHS_mc.root
  TString rootFile[2]; 
  rootFile[0]  = pathName + (TString) "hPhoton1Pt_PFCHS_mc.root";
  rootFile[1]  = pathNameQCD + (TString) "hPhoton1Pt_PFCHS_data.root"; 

  TFile* file[2];
  TH1D* histo[2];
  TH1D* onlyGamma;
  TH1D* ratio;
  
  for(int i =0; i<2; i++){
    file[i] = TFile::Open(rootFile[i]);
    file[i]->GetObject("histo",histo[i]); 
    histo[i] ->Rebin(20);   
  }

  file[0] = TFile::Open(rootFile[0]);
  file[0]->GetObject("histo",onlyGamma); 
  onlyGamma ->Rebin(20);

  //histo[1]->Add(histo[0],-1);

  onlyGamma->Add(histo[1],1);


  TCanvas *stackPlot = new TCanvas("stackPlot","stackPlot",500,500,500,500);

  TLegend *l = new TLegend(0.5,0.8,0.9,0.9);


 
 

  stackPlot->SetLogy();
  onlyGamma->GetXaxis()->SetRangeUser(0.,800.);
  histo[1]->GetXaxis()->SetRangeUser(0.,800.);
  onlyGamma->SetMaximum(pow(10.,8));
  onlyGamma->SetMinimum(0.1);
  histo[1]->SetMaximum(pow(10.,8));
  histo[1]->SetMinimum(0.1);
  onlyGamma->SetFillColor(2);
  histo[1]->SetFillColor(5);
  onlyGamma->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  onlyGamma->Draw("hist");
  histo[1]->Draw("samehist");

 l -> SetFillColor(0);
  l -> SetTextSize(0.042);
  l -> AddEntry(onlyGamma,"#gamma + jet MC","f");
  l -> AddEntry(histo[1],"QCD","f");

  l->Draw("SAME");

  stackPlot->SaveAs("plotsQCD/stackGammaJetQCD.pdf");
  
  file[1] = TFile::Open(rootFile[1]);
  file[1]->GetObject("histo",ratio); 
  ratio->Rebin(20);  
   
  double correlation = 0;

  correlation = sqrt(histo[0]->GetSumOfWeights()/histo[1]->GetSumOfWeights());

  cout<<"correlation is = "<<correlation<<endl<<endl<<endl;

  cout<<"wo QCD  : histo[0]->GetSumOfWeights() = "<<histo[0]->GetSumOfWeights()<<endl;
  cout<<"wo QCD  : histo[0]->Integral()        = "<<histo[0]->Integral()<<endl;
  cout<<"wo QCD  : histo[0]->GetEntries()      = "<<histo[0]->GetEntries()<<endl;
  cout<<"with QCD: histo[1]->GetSumOfWeights() = "<<histo[1]->GetSumOfWeights()<<endl;
  cout<<"with QCD: histo[1]->Integral()        = "<<histo[1]->Integral()<<endl;
  cout<<"with QCD: histo[1]->GetEntries()      = "<<histo[1]->GetEntries()<<endl;
  cout<<"only QCD:         Integral            = "<<histo[1]->Integral()-histo[0]->Integral()<<endl;
  cout<<"only QCD:         Entries             = "<<histo[1]->GetEntries()-histo[0]->GetEntries()<<endl<<endl;

  ratio -> Add(histo[0],-1);
  ratio -> Divide(histo[1]);

  gStyle -> SetPadTopMargin(0.15);
  TCanvas *plotRatio = new TCanvas("ratio","ratio",500,500,500,500);

  plotRatio->cd();
  ratio->SetMinimum(0.0);
  ratio->GetXaxis()->SetRangeUser(40.,800.);
  ratio -> SetTitle("#int(QCD)/#int(QCD +(#gamma+jet))");
  ratio->Draw();
  tot_filename = (TString) "plotsQCD/qcdContamination.pdf";
  plotRatio -> SaveAs(tot_filename);
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
  rootFile[0]  = pathName + (TString) "Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_mc_" +method + (TString) ".root";
  rootFile[1]  = pathNameQCD + (TString) "Resolution_for_" + (long) eta + (TString) "_eta_bin_PFCHS_data_" +method + (TString)".root";

  etaString = "Uncertainty on QCD contamination";
  
  mg -> SetTitle(etaString);
    

  TCanvas *c = new TCanvas("c",etaString,200,10,800,800);
  c -> cd();

  TString fitName;

  double maximum = 0;

  TGraphErrors* graph[2]; 
  for(int i =0; i<2; i++){
    graph[i] = readTGraphErrors(rootFile[i],"Graph","Graph");    
    graph[i] -> SetMarkerStyle(20);
    if(graph[i]->GetYaxis()->GetXmax() > maximum) maximum = graph[i]->GetYaxis()->GetXmax();
  }
  
  graph[0] -> SetMarkerColor(8);
  graph[1] -> SetMarkerColor(9);
 
  graph[0] -> SetLineColor(8);
  graph[1] -> SetLineColor(9);

  graph[0] -> GetFunction("fResolution") -> SetLineColor(8);
  graph[1] -> GetFunction("fResolution") -> SetLineColor(9);
  
 
  TLegend *legend  = new TLegend(0.5,0.8,0.9,0.9);
  legend -> SetFillColor(0);
  legend -> SetTextSize(0.033);
  legend -> AddEntry(graph[0],"without QCD","l");
  legend -> AddEntry(graph[1],"with QCD","l");
  
  mg->Add(graph[0]);  
  mg->Add(graph[1]);  
  
  mg->Draw("AP");
    
  mg -> GetYaxis() -> SetTitle("JER");
  mg -> GetYaxis() -> SetTitleOffset(1.3); 
  mg -> GetYaxis() -> SetTitleSize(0.05);
  mg -> SetMinimum(0.02);
  mg -> SetMaximum(0.18);   
  mg -> GetXaxis() -> SetLimits(0,600);

  mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  mg -> GetXaxis() -> SetTitleSize(0.05);
  

  TLatex*  info   = new TLatex();
  info->SetTextFont(132);
  info-> SetNDC();
  info->SetTextSize(0.045);
  AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBinsQCD[eta-1],etaBinsQCD[eta]);
  info->DrawLatex(0.6,0.7,AuxString);

  
  legend->Draw("same");
  tot_filename = (TString) "plotsQCD/Resolution_for_" + (long) eta + (TString) "_eta_bin_QCDUncertainty_" + method + (TString) ".pdf";
  c -> SaveAs(tot_filename);

  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // 2.) Relative uncertainty to MC -> MC*(1 +- Delta)
  double *dataY, *dataQCDY, *dataX, *dataQCDX, *dataEY, *dataQCDEY, *dataEX, *dataQCDEX;

  double *errorQCDY = new double[graph[0]->GetN()];
  double *errorQCDEY= new double[graph[0]->GetN()];

  dataY     = graph[0] -> GetY();
  dataQCDY  = graph[1] -> GetY();
  dataX     = graph[0] -> GetX();
  dataQCDX  = graph[1] -> GetX();
  dataEY    = graph[0] -> GetEY();
  dataQCDEY = graph[1] -> GetEY();
  dataEX    = graph[0] -> GetEX();
  dataQCDEX = graph[1] -> GetEX();

  int nData = graph[0]->GetN();

  if(graph[0]->GetN() != graph[1]->GetN()) cout<<endl<<"The number of entries in the various graphs are not the same (diff = "<<abs(graph[1]->GetN()-graph[0]->GetN())<<")!!!!!!!!!!!!!!!!  (in "<<rootFile[1]<<") "<<endl<<endl;

  //cout<<endl<<"Entries in "<<rootFile[0]<<" = "<<graph[0]->GetN()<<endl;
  //cout<<"Lowest Pt Bin = "<<dataX[0]<<endl;
  //cout<<"Highest Pt Bin = "<<dataX[graph[0]->GetN()-1]<<endl<<endl;

  //cout<<"Entries in "<<rootFile[1]<<" = "<<graph[1]->GetN()<<endl;
  //cout<<"Lowest Pt Bin (low) = "<<dataQCDX[0]<<endl;
  //cout<<"Highest Pt Bin = "<<dataQCDX[graph[1]->GetN()-1]<<endl<<endl;


  for(int i=0; i<nData; i++){

    errorQCDY[i]  = dataQCDY[i]/dataY[i] - 1;
    errorQCDEY[i] = TMath::Sqrt(TMath::Power(1./dataY[i],2)*TMath::Power(dataQCDEY[i],2) + TMath::Power(dataQCDY[i]/TMath::Power(dataY[i],2),2)*TMath::Power(dataEY[i],2) - correlation*2.0 * dataQCDY[i]/TMath::Power(dataY[i],3)*dataQCDEY[i]*dataEY[i]);

    dataX[i] = dataX[i];
    dataEX[i] = dataEX[i];
  }

  TGraphErrors *plotQCD = new TGraphErrors(nData,dataX,errorQCDY,dataEX,errorQCDEY);

  // Search for the interval where 68% of all point are included (= 11 points)
  double* arraySort = new double[nData];
  float auxSort = 10000000.;

  
  for(int i=0; i<nData; i++){
    arraySort[i]   = abs(errorQCDY[i]);
  }

  for(int j=0; j<nData; j++){

    for(int i=0; i<nData-1; i++){

      if(arraySort[i] > arraySort[i+1]){
	auxSort        = arraySort[i];
	arraySort[i]   = arraySort[i+1];
	arraySort[i+1] = auxSort;
      }
    }
  }

  TF1* sigma1Interval = new TF1("sigma1Interval","pol0",0,600);
  sigma1Interval->SetLineColor(8);
  
  if(2*arraySort[5]<arraySort[7]){
    sigma1Interval->SetParameter(0,arraySort[7]/2.);
    cout<<"arraySort[5]    = "<<arraySort[5]<<endl;
    cout<<"arraySort[7]    = "<<arraySort[7]<<endl;
    cout<<"arraySort[7]/2. = "<<arraySort[7]/2.<<endl;
    cout<<"arraySort[8]    = "<<arraySort[8]<<endl;
  }
  else  sigma1Interval->SetParameter(0,arraySort[5]);

  TF1* sigma1Intervaldown = new TF1("sigma1Intervaldown","pol0",0,600);
  sigma1Intervaldown->SetParameter(0,-sigma1Interval->GetParameter(0));
  sigma1Intervaldown->SetLineColor(8);
  

  
  TCanvas *c1 = new TCanvas("c1","c1",200,10,800,800);
  c1 -> cd();
  c1 -> SetBottomMargin(0.14);
  c1 -> SetLeftMargin(0.14);

  plotQCD -> SetMarkerStyle(24);  
  plotQCD -> SetMarkerColor(9);  
  plotQCD -> SetLineColor(9); 
 
  TLegend *legend1  = new TLegend(0.5,0.8,0.9,0.9);
  legend1 -> SetFillColor(0);
  legend1 -> SetTextSize(0.033);
  legend1 -> AddEntry(plotQCD,"with QCD","p");

  mg = new TMultiGraph();
  mg->Add(plotQCD);  
 

  etaString = "Relative QCD Uncertainty";

  mg -> SetTitle(etaString);  
  
  mg->Draw("AP");
  sigma1Interval->Draw("same");
  sigma1Intervaldown->Draw("same");      

  mg -> GetYaxis() -> SetTitle("JER_{with QCD} /JER_{without QCD} ");
  mg -> GetYaxis() -> SetTitleOffset(1.3); 
  mg -> GetYaxis() -> SetTitleSize(0.05);
  mg -> SetMinimum(-0.50);
  mg -> SetMaximum(0.50); 
  mg -> GetXaxis() -> SetLimits(0,600);  
  mg -> GetXaxis() -> SetTitleSize(0.05);
  
  mg -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  //mg -> GetXaxis() -> SetTitleOffset(1.2); 
  
  legend1->Draw("same");

  TLatex*  info1   = new TLatex();
  info1->SetTextFont(132);
  info1-> SetNDC();
  info1->SetTextSize(0.045);
  AuxString = Form("%4.1f < |#eta^{Jet}| < %4.1f",etaBinsQCD[eta-1],etaBinsQCD[eta]);
  info1->DrawLatex(0.6,0.7,AuxString);

  TF1 *line = new TF1("line","pol0",0,600);
  line -> SetLineColor(1);
  line ->SetParameter(0,0.);
  line->Draw("same");

  TF1 *fitLine = new TF1("fitLine","pol0",0,600);
  //TF1 *fitLine = new TF1("fitLine","expo",0,600);
  fitLine -> SetLineColor(8);
  plotQCD -> Fit("fitLine","QR");

  TString text = "#scale[0.9]{68% of points are included in Interval:} " + (TString) Form("#scale[0.9]{f = %4.3f}",sigma1Interval->GetParameter(0)) ;
  info1->DrawLatex(0.2,0.3,text);

  AuxString = Form("Chi2/ndof = %4.1f/%i",fitLine->GetChisquare(),fitLine->GetNDF());
  info1->DrawLatex(0.5,0.2,AuxString);

  AuxString = Form("f = %4.3f #pm %4.3f",fitLine->GetParameter(0),fitLine->GetParError(0));
  info1->DrawLatex(0.5,0.25,AuxString);

  tot_filename = (TString) "plotsQCD/Relative_Resolution_for_" + (long) eta + (TString) "_eta_bin_QCDUncertainty_" + method + (TString) ".pdf";
  c1 -> SaveAs(tot_filename);

  
  TF1* fResolutionWoQCD   = graph[0]->GetFunction("fResolution");
  TF1* fResolutionWithQCD = graph[1]->GetFunction("fResolution");
  //fResolutionWoQCD   -> SetName("fResolutionWoQCD");
  //fResolutionWithQCD -> SetName("fResolutionWithQCD");


  TF1* diff = new TF1("diff","TMath::Sqrt(TMath::Sign(1.,[0])*TMath::Power([0]/x,2)+TMath::Power([1],2)*TMath::Power(x,[3]-1)+TMath::Power([2],2))/(TMath::Sqrt(TMath::Sign(1.,[4])*TMath::Power([4]/x,2)+TMath::Power([5],2)*TMath::Power(x,[7]-1)+TMath::Power([6],2))) - 1.",40,600);

  diff->SetParameter(0,fResolutionWithQCD->GetParameter(0));
  diff->SetParameter(1,fResolutionWithQCD->GetParameter(1));
  diff->SetParameter(2,fResolutionWithQCD->GetParameter(2));
  diff->SetParameter(3,fResolutionWithQCD->GetParameter(3));
  diff->SetParameter(4,fResolutionWoQCD->GetParameter(0));
  diff->SetParameter(5,fResolutionWoQCD->GetParameter(1));
  diff->SetParameter(6,fResolutionWoQCD->GetParameter(2));
  diff->SetParameter(7,fResolutionWoQCD->GetParameter(3));


  TCanvas *c4 = new TCanvas("c4","c4",200,10,800,800);
  c4 -> cd();
  c4 -> SetBottomMargin(0.14);
  c4 -> SetLeftMargin(0.14);

  
  
  diff -> Draw();

  tot_filename = (TString) "plotsQCD/FinalErrorsQCD_" + type + (TString) "_" + method + (TString) ".root";
  TFile *f = new TFile(tot_filename,"RECREATE");
  f -> WriteTObject(diff,"function");
  f->Close();
  delete f;

  diff -> GetYaxis() -> SetTitle("#Delta JER /JER ");
  diff -> GetYaxis() -> SetTitleOffset(1.3); 
  diff -> SetMinimum(-0.05);
  diff -> SetMaximum(0.2); 
  diff -> GetXaxis() -> SetLimits(0,600);  
  diff -> GetXaxis() -> SetTitleSize(0.05);
  diff -> GetYaxis() -> SetTitleSize(0.05);  
  diff -> GetXaxis() -> SetTitle("p_{T}^{#gamma} [GeV]");
  diff -> SetTitle("");

  tot_filename = (TString) "plotsQCD/FinalErrorsQCD_" + type + (TString) "_" + method + (TString) ".pdf";
  c4->SaveAs(tot_filename);

  return 0;

}
