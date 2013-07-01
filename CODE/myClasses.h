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
  double alpha[2*alpha_int], alphaError[2*alpha_int], mean[2*alpha_int], meanError[2*alpha_int], sigma[2*alpha_int], sigmaError[2*alpha_int];
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
  double pT[pt_int], pTError[pt_int], scale[pt_int], scaleError[pt_int], resolution[pt_int], resolutionError[pt_int], q[pt_int], qError[pt_int], qprime[pt_int], qprimeError[pt_int];
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

CResponse* JetResponsePhotonHemisphere[pt_int][eta_int][alpha_int];
CResponse* JetResponseJetHemisphere[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsic[pt_int][eta_int][alpha_int];
CResponse* JetImbalancePhotonHemisphere[pt_int][eta_int][alpha_int];
CResponse* JetImbalanceJetHemisphere[pt_int][eta_int][alpha_int];
CResponse* JetResponse[pt_int][eta_int][alpha_int];

CResponse* JetIntrinsicle5[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsicle10gt5[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsicle15gt10[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsicle20gt15[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsicgt20[pt_int][eta_int][alpha_int];

CScaleResAlpha* JetScaleResAlpha[pt_int][eta_int];
CScaleResAlpha* JetIntrinsicAlpha[pt_int][eta_int];
CScaleResAlpha* JetImbalanceAlpha[pt_int][eta_int];


CScaleRes* JetScaleRes[eta_int];
CScaleRes* JetImbalancePt[eta_int];
CScaleRes* JetIntrinsicPt[eta_int];

CResponse* JetResponseGluon[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsicGluon[pt_int][eta_int][alpha_int];
CResponse* JetImbalanceGluon[pt_int][eta_int][alpha_int];

CResponse* JetResponseQuark[pt_int][eta_int][alpha_int];
CResponse* JetIntrinsicQuark[pt_int][eta_int][alpha_int];
CResponse* JetImbalanceQuark[pt_int][eta_int][alpha_int];

// For PU uncertainty
CResponse* JetResponsePUle10[pt_int][eta_int][alpha_int];
CResponse* JetResponsePUgt10le15[pt_int][eta_int][alpha_int];
CResponse* JetResponsePUgt15le20[pt_int][eta_int][alpha_int];
CResponse* JetResponsePUgt20[pt_int][eta_int][alpha_int];


CScaleResAlpha* JetScaleResAlphaGluon[pt_int][eta_int];
CScaleResAlpha* JetIntrinsicAlphaGluon[pt_int][eta_int];
CScaleResAlpha* JetImbalanceAlphaGluon[pt_int][eta_int];

CScaleResAlpha* JetScaleResAlphaQuark[pt_int][eta_int];
CScaleResAlpha* JetIntrinsicAlphaQuark[pt_int][eta_int];
CScaleResAlpha* JetImbalanceAlphaQuark[pt_int][eta_int];

CScaleRes* JetScaleResGluon[eta_int];
CScaleRes* JetImbalancePtGluon[eta_int];
CScaleRes* JetIntrinsicPtGluon[eta_int];

CScaleRes* JetScaleResQuark[eta_int];
CScaleRes* JetImbalancePtQuark[eta_int];
CScaleRes* JetIntrinsicPtQuark[eta_int];

//CsysUncerMC* sysUncer[eta_int];

JetIndexCol corrJets;

#endif
