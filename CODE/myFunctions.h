#ifndef myFunctions_h
#define myFunctions_h


#include "TString.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"


void plotTH1(TH1* histo,const TString &filename, bool logy);

void plotTH2(TH2* histo,const TString &filename);

void plotTGraphErrors(TGraphErrors* graph, TString &filename, TString &graphTitle, TString &xTitle, TString &yTitle, double x_low, double x_up, double y_low, double y_up, char legEntry[100], bool legend);

void saveObject(TObject* object, TString filename, TString objectName);


TH1D* createTH1(TH1D *histo, const TString &HistoName,const TString &HistoTitle,int nBin, double bound_low, double bound_up,const TString &xTitle);

TH2D*createTH2(TH2D *histo,const TString &HistoTitle,int XnBin,double Xbound_low,double Xbound_up,int YnBin,double Ybound_low,double Ybound_up,const TString &xTitle,const TString &yTitle);

TH1D* readTH1(const TString &fileName, const TString &histName, const TString &newHistName);

TGraphErrors* readTGraphErrors(const TString &fileName, const TString &gName, const TString &newGName);

#endif
