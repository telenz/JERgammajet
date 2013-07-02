#ifndef myClasses_h
#define myClasses_h

#include "myDeclarations.h"
#include <iostream>

#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"

//----
class CResponse{

public:
  
  TH1D *hResponse, *hAlpha, *hPt, *hPtLeadingJet, *hgenPtLeadingJet;
  double sigma_array, sigmaAUX, sigma_error, mean_array, mean_error, alpha_array, alpha_error;
  double pT_array, pT_error, RMS_array, RMS_error, chi2, chi2_ndof, meanAUX, meanAUX2;
  int ndof;
  bool imbalance;
  float sigma;
  double bin[2];
  double prob[2];
  
    
  CResponse(int typeResponse);
 
  ~CResponse();

  void calculatePt();

  void calculate();
  };

//----

class CScaleResAlpha{

public: 
  TF1 *fScaleAlpha, *fResolutionAlpha;
  TGraphErrors *gJetScaleAlpha, *gJetResolutionAlpha ;
  double pT, pTError, scale, scaleError, resolution, resolutionError;   
  double alpha[2*nAlphaBins], alphaError[2*nAlphaBins], mean[2*nAlphaBins], meanError[2*nAlphaBins], sigma[2*nAlphaBins], sigmaError[2*nAlphaBins];
  float q, qprime, c, cprime, cprimeError, m, mprime, mprimeError, qError, qprimeError;
    
  CScaleResAlpha();

  void calculate(int length, int fit);
    
  ~CScaleResAlpha();

};
//----
class CScaleRes{

   
  public:
  TF1* fScale, *fResolution;
  TGraphErrors* gScale, *gResolution, *gq, *gqprime;
  double scalenumber, scalenumberError, N, NErr, S, SErr, C, CErr, m, mErr;
  double pT[nPtBins], pTError[nPtBins], scale[nPtBins], scaleError[nPtBins], resolution[nPtBins], resolutionError[nPtBins], q[nPtBins], qError[nPtBins], qprime[nPtBins], qprimeError[nPtBins];
  double maximum, minimum;
    
  CScaleRes();

  void calculate(int length);
 
  ~CScaleRes();
};

class JetIndexCol {
private:
  class Jet {
  public:
    Jet(unsigned int jetIdx, double jetMomentum) : idx_(jetIdx), pt_(jetMomentum) {};
    const unsigned int idx_;
    const double pt_;
    // For sorting jets in pt
    static bool ptGreaterThan(const Jet *idx1, const Jet *idx2) {
      // check for 0
      if(idx1 == 0) {
	return idx2 != 0;
      } else if (idx2 == 0) {
	return false;
      } else {
	    return idx1->pt_ > idx2->pt_;
      }
    }
  };
  
  std::vector<Jet*> jets_;
  
  
public:
  JetIndexCol() {}
  ~JetIndexCol() { clear(); }
  
  unsigned int operator()(unsigned int i) { return idx(i); }
  unsigned int nJets() const { return jets_.size(); }
  unsigned int idx(unsigned int i) const { return jets_.at(i)->idx_; }
  double pt(unsigned int i) const { return jets_.at(i)->pt_; }
  
  void add(unsigned int jetIdx, double jetMomentum) {
    jets_.push_back(new Jet(jetIdx,jetMomentum));
  }
  void clear() {
    for(std::vector<Jet*>::iterator it = jets_.begin(); it != jets_.end(); ++it) {
      //cout<<"it = "<<*it<<endl;
      delete *it;
      //cout<<"after deleting"<<endl;
    }
    jets_.clear();
  }
  void sort() {
    std::sort(jets_.begin(),jets_.end(),Jet::ptGreaterThan);
  }
};

// Declarations of class objects

CResponse* JetResponsePhotonHemisphere[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetResponseJetHemisphere[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsic[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetImbalancePhotonHemisphere[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetImbalanceJetHemisphere[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetResponse[nPtBins][nEtaBins][nAlphaBins];

CResponse* JetIntrinsicle5[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsicle10gt5[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsicle15gt10[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsicle20gt15[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsicgt20[nPtBins][nEtaBins][nAlphaBins];

CScaleResAlpha* JetScaleResAlpha[nPtBins][nEtaBins];
CScaleResAlpha* JetIntrinsicAlpha[nPtBins][nEtaBins];
CScaleResAlpha* JetImbalanceAlpha[nPtBins][nEtaBins];


CScaleRes* JetScaleRes[nEtaBins];
CScaleRes* JetImbalancePt[nEtaBins];
CScaleRes* JetIntrinsicPt[nEtaBins];

CResponse* JetResponseGluon[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsicGluon[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetImbalanceGluon[nPtBins][nEtaBins][nAlphaBins];

CResponse* JetResponseQuark[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetIntrinsicQuark[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetImbalanceQuark[nPtBins][nEtaBins][nAlphaBins];

// For PU uncertainty
CResponse* JetResponsePUle10[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetResponsePUgt10le15[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetResponsePUgt15le20[nPtBins][nEtaBins][nAlphaBins];
CResponse* JetResponsePUgt20[nPtBins][nEtaBins][nAlphaBins];


CScaleResAlpha* JetScaleResAlphaGluon[nPtBins][nEtaBins];
CScaleResAlpha* JetIntrinsicAlphaGluon[nPtBins][nEtaBins];
CScaleResAlpha* JetImbalanceAlphaGluon[nPtBins][nEtaBins];

CScaleResAlpha* JetScaleResAlphaQuark[nPtBins][nEtaBins];
CScaleResAlpha* JetIntrinsicAlphaQuark[nPtBins][nEtaBins];
CScaleResAlpha* JetImbalanceAlphaQuark[nPtBins][nEtaBins];

CScaleRes* JetScaleResGluon[nEtaBins];
CScaleRes* JetImbalancePtGluon[nEtaBins];
CScaleRes* JetIntrinsicPtGluon[nEtaBins];

CScaleRes* JetScaleResQuark[nEtaBins];
CScaleRes* JetImbalancePtQuark[nEtaBins];
CScaleRes* JetIntrinsicPtQuark[nEtaBins];

JetIndexCol corrJets;

#endif
