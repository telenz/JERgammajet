#include "myDeclarations.h"
#include "myClasses.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"

bool applyCuts(){


  // Check if array is large enough
  if( nobjJet > NJETS || nobjPhoton > NPHOTONS) {
    std::cerr << "\nERROR: 'nobjJet = " << nobjJet << " > NJETS = " << NJETS << "'" << std::endl;
    std::cerr << "\nERROR: 'nobjPhoton = " << nobjJet << " > NPHOTONS = " << NPHOTONS << "'" << std::endl;
    cut0 += 1;
    return 0;
  }

  idx1stJet = 0;
  idx2ndJet = 0;
  alpha     = 0;

  //------------------------------------------------------------------------------(1. CUT)------------------------
  // Select events with at least 1 Photon         (1. CUT) 
  if( nobjPhoton < 1 ){
    cut1 = cut1 +1;
    //cout<<"Discarded because there is no photon in this event!"<<endl;
    return 0;}
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(2. CUT)------------------------
  // Select only events with tight leading Photon (2. CUT)
  if(!tight[0]){
    cut2 = cut2 +1;
    //cout<<"Discarded because leading Photon is not tight!"<<endl;
    return 0;
  }
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(3. CUT)------------------------
  // Photon Has Pixel Seed                       (3. CUT)
  if(photonHasPixelSeed[0]){
    cut3 = cut3 +1;
    //cout<<"Discarded because of PhotonHasPixelSeed!"<<endl;
    return 0;
  }
  //--------------------------------------------------------------------------------------------------------------

  //**************************************************************************************************************
  // Correct Jets and Order them 
  corrJets.clear();
  for(int j = 0; j < nobjJet; j++) corrJets.add(j,/*(1.+jetUncert[j])*/jetCorrL1[j]*jetCorrL2L3[j]*jetPt[j]);
  corrJets.sort();
  
  // Definition of some variables
  double_t deta[nobjJet];
  double_t dR[nobjJet];
  double_t dphi[nobjJet];
  
  int j_jet    = 0;
  int j_photon = 0;
  idx1stJet     = -1;
  idx2ndJet        = -1;
  idxPhoton    = -1;
    
  // Find out to what index the photon belongs to in the jet sample
  for (int i=0 ; i<nobjJet ; i++) {
    
    if(j_jet==2 && j_photon==1) continue;

    dphi[i] = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(i)]-photonPhi[0]));
    deta[i] = std::abs(jetEta[corrJets.idx(i)] - photonEta[0]);
    dR[i]   = TMath::Sqrt(dphi[i]*dphi[i] + deta[i]*deta[i]);

    if(dR[i] > 0.5 && j_jet<2) {
      
      if(j_jet == 0){
        idx1stJet = i;
        j_jet    = j_jet+1;
      }       
      else if(j_jet == 1){
	idx2ndJet = i;
	j_jet = j_jet+1;
      }
      else cout<<"something wrong here"<<endl;      
    }
    else if(dR[i]<=0.5){
     
      idxPhoton = i; 
      j_photon  = j_photon+1;              
    }
  }
  //**************************************************************************************************************

  //------------------------------------------------------------------------------(4. CUT)------------------------
  // Discard events with no leading Jet in the sample   (4. CUT)
  if(idx1stJet==-1) {
    cut4 = cut4 +1;
    //cout<<"Discarded because there is no balancing Jet in the Jet Sample!"<<endl;
    return 0;   
  }
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(5. CUT)------------------------
  // Take only tight Jets                               (5. CUT)
  if(!jetIDTight[corrJets.idx(idx1stJet)]){
    cut5 = cut5 +1;
    // cout<<"Discarded because of tight Jet ID "<<endl;
    return 0;}
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(6. CUT)------------------------
  // Reject events with |JetEta| > 5 or jetPt < 10      (6. CUT)
  if(std::abs(jetEta[corrJets.idx(idx1stJet)])>5.0 || jetPt[corrJets.idx(idx1stJet)] < 10.){  
    cut6 = cut6 +1;
    //cout<<"Discarded because of Jet eta!"<<endl; 
    return 0;
  }
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(7. CUT)------------------------
  // Upper and lower Bounds on all Binning variables    (7. CUT)
  if(photonPt[0] < ptBins[0] || photonPt[0] > ptBins[nPtBins] || std::abs(jetEta[corrJets.idx(idx1stJet)]) > etaBins[nEtaBins]){
    cut7 = cut7 +1;
    return 0;
  }
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(8. CUT)-----------------------
  // photonEta > 1.3 Cut                                (8. CUT)
  if(std::abs(photonEta[0])>1.3){  
    cut8 = cut8 +1;  
    //cout<<"Discarded because of photon etacut!"<<endl; 
    return 0;}
  //--------------------------------------------------------------------------------------------------------------
  
  //Get rid of strange events (just for 2011 data)
  if(date == 2011 && idxPhoton!=-1 && corrJets.pt(idxPhoton)/photonPt[0]<0.5) return 0;

  //**************************************************************************************************************
  // Find the respective generator Jet 
  gen1stJetidx = -20000000;
  gen2ndJetidx = -20000000;

  jetPt1stJet  = log(0.);

  float cSmearing = 1.;
  int   test      = 0;
  bool  jet1      = true;
  bool  jet2      = true;
  
  if(isMC || testClosure){
    for (int i=0 ; i<nobjGenJet ; i++) {
      
      if(genJetColJetIdx[i] == corrJets.idx(idx1stJet) && jet1){
      	gen1stJetidx = i;
	test += 1;
	jet1 = false;
      }
      
      if(idx2ndJet != -1){
	if(genJetColJetIdx[i] == corrJets.idx(idx2ndJet) && jet2){
	  gen2ndJetidx = i;
	  test += 1;
	  jet2 = false;
	}
      }
      
      if(test == 2) break;      
    }
  
  } 
  //**************************************************************************************************************
  
    
  //------------------------------------------------------    
  if(isMC || testClosure){

    jetPt1stJet = corrJets.pt(idx1stJet);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // MC smearing Code piece
    //cout<<"Smearing is activated!!"<<endl;
    
    // Smear Second Jet
    if(testClosure && !isMC){
      //cout<<"You are Smearing all Jets!!"<<endl;
      if(idx2ndJet != -1 && gen2ndJetidx >= 0){
	if(std::abs(jetEta[corrJets.idx(idx2ndJet)])<0.5) cSmearing = 1.05;
	else if(std::abs(jetEta[corrJets.idx(idx2ndJet)])<1.1 && std::abs(jetEta[corrJets.idx(idx2ndJet)])>0.5) cSmearing = 1.07;
	else if(std::abs(jetEta[corrJets.idx(idx2ndJet)])<1.7 && std::abs(jetEta[corrJets.idx(idx2ndJet)])>1.1) cSmearing = 1.09;
	else if(std::abs(jetEta[corrJets.idx(idx2ndJet)])<2.3 && std::abs(jetEta[corrJets.idx(idx2ndJet)])>1.7) cSmearing = 1.11;
	else cSmearing = 1.000;
	jetPt2nd = genJetPt[gen2ndJetidx] + cSmearing*(jetPt2nd - genJetPt[gen2ndJetidx]);
      }

      // Smear first Jet
      if(std::abs(jetEta[corrJets.idx(idx1stJet)])<0.5) cSmearing = 1.05;
      else if(std::abs(jetEta[corrJets.idx(idx1stJet)])<1.1 && std::abs(jetEta[corrJets.idx(idx1stJet)])>0.5) cSmearing = 1.07;
      else if(std::abs(jetEta[corrJets.idx(idx1stJet)])<1.7 && std::abs(jetEta[corrJets.idx(idx1stJet)])>1.1) cSmearing = 1.09;
      else if(std::abs(jetEta[corrJets.idx(idx1stJet)])<2.3 && std::abs(jetEta[corrJets.idx(idx1stJet)])>1.7) cSmearing = 1.11;
      else cSmearing = 1.000;
      
      if(gen1stJetidx >= 0) jetPt1stJet = genJetPt[gen1stJetidx] + cSmearing*(corrJets.pt(idx1stJet) - genJetPt[gen1stJetidx]);
      else return 0;
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(gen1stJetidx >= 0) intrinsic   = jetPt1stJet/genJetPt[gen1stJetidx];   
    if(gen1stJetidx >= 0) imbalance   = genJetPt[gen1stJetidx]/photonPt[0]; 
  }
  else if(!isMC) jetPt1stJet = corrJets.pt(idx1stJet);

  
  //------------------------------------------------------------------------------(9. CUT)-----------------------
  // Cut on second Jetpt > 10 GeV                   (9. CUT)
  if(idx2ndJet != -1){
    
    if(corrJets.pt(idx2ndJet) > 10.) jetPt2nd = corrJets.pt(idx2ndJet);
    else{
      cut9 += 1;
      return 0;
    }
  }
  else{
    
    jetPt2nd = 0.;
    no2ndJetinEvent += 1;
  }
  //-------------------------------------------------------------------------------------------------------------

  // Calculate response    
  response       = jetPt1stJet/photonPt[0]; 
  float deltaphi = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets.idx(idx1stJet)]-photonPhi[0]));
 
  /*
    if(idx2ndJet != -1){
    photonVector.SetPtEtaPhiE(photonPt[0],photonEta[0],photonPhi[0],photonE[0]);
    jet2ndVector.SetPtEtaPhiE(jetPt[corrJets.idx(idx2ndJet)],jetEta[corrJets.idx(idx2ndJet)],jetPhi[corrJets.idx(idx2ndJet)],jetE[corrJets.idx(idx2ndJet)]); 
    sumVector = photonVector + jet2ndVector;
    }
  
  
    for(int i=0; i<nPtBins-1; i++){ 
	      
    if(photonPt[0] >= ptBins[i] && photonPt[0] < ptBins[i+1]){
    alpha = jetPt2nd/((ptBins[i]+ptBins[i+1])/2.) * 100.;  
    }
    else if(photonPt[0]>ptBins[nPtBins-1]) alpha = jetPt2nd/(500.) * 100.;
    }
  */

  //alpha = jetPt2nd/((jetPt1stJet+photonPt[0])/2.) * 100.;
  //alpha = jetPt2nd/jetPt1stJet * 100.;
  //alpha = jetPt2nd/(sumVector.Pt()) * 100.;
  alpha = jetPt2nd/photonPt[0] * 100.;
      
  //------------------------------------------------------------------------------(10. CUT)-----------------------
  // Cut on highest Alpha                             (10. CUT)
  if(idx2ndJet!=-1){
    if(jetPt2nd > alphaBins[nAlphaBins]/100. * photonPt[0]){
      cut10=cut10+1;
      return 0;
      // cout<<"Discarded because of alpha constraint. alpha = "<<alpha <<endl;
    }
  }
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(11. CUT)------------------------
  // Select events with back-to-back photon and jet   (11. CUT)
  if( deltaphi < deltaPhiCutValue){
    cut11 = cut11 +1;
    //cout<<"Discarded because of delta phi!"<<endl;
    return 0;
  }
  //--------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------(12. CUT)-----------------------
  // Cuts on Trigger                                  (12. CUT)
  if((!isMC && !QCDUncertaintyEvaluation) || applyTriggeronMC ){

    if((photonPt[0] >= ptBins[0] && photonPt[0] < ptBins[1]) && !hltPhoton[0]){
      cut12[0] += 1;
      return 0;}
    if((photonPt[0] >= ptBins[1] && photonPt[0] < ptBins[2]) && !hltPhoton[1]){
      cut12[1] += 1; 
      return 0;}
    if((photonPt[0] >= ptBins[2] && photonPt[0] < ptBins[3]) && !hltPhoton[2]){
      cut12[2] += 1; 
      return 0;}
    if((photonPt[0] >= ptBins[3] && photonPt[0] < ptBins[4]) && !hltPhoton[3]){
      cut12[3] += 1;
      return 0;}
    if((photonPt[0] >= ptBins[4] && photonPt[0] < ptBins[5]) && !hltPhoton[4]){
      cut12[4] += 1;
      return 0;}
    if(date == 2012){
      if((photonPt[0] >= ptBins[5] && photonPt[0] < ptBins[6]) && !hltPhoton[5]){
	cut12[5] += 1; 
	return 0;}
      if((photonPt[0] >= ptBins[6]) && !hltPhoton[6]){
	cut12[6] += 1;
	return 0;}
    }
  }  
  //--------------------------------------------------------------------------------------------------------------
 
  //------------------------------------------------------------------------------(13. CUT)------------------------
  // Only for Flavor Unceratinty (13. CUT)
  bool flavorFullfilled = true;
  if(isMC && flavorUncertaintyEvaluationPhys){
    if(setFlavorSelection == 1){      if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 0) flavorFullfilled = false;}
    else if(setFlavorSelection == 2){ if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) > 5 || abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 0) flavorFullfilled = false;}
    else if(setFlavorSelection == 3){ if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) != 21) flavorFullfilled = false;}
    if(!flavorFullfilled){
      cut14 = cut14 +1;
      return 0;
    }
  }
  if(isMC && flavorUncertaintyEvaluationAlgo){
    if(setFlavorSelection == 1){      if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 0) flavorFullfilled = false;}
    else if(setFlavorSelection == 2){ if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) > 5 || abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 0) flavorFullfilled = false;}
    else if(setFlavorSelection == 3){ if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) != 21) flavorFullfilled = false;}
    else if(setFlavorSelection == 4){ if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) > 3 || abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 0) flavorFullfilled = false;}
    else if(setFlavorSelection == 5){ if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) > 5 || abs(genJetID_algo[corrJets.idx(idx1stJet)]) < 4)  flavorFullfilled = false;}
    if(!flavorFullfilled){
      cut14 = cut14 +1;
      return 0;
    }
  }
  //--------------------------------------------------------------------------------------------------------------

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Make plots after whole selection
 
  // Number of Vertices for different pt bins  
  for(int i=0; i<numTrigger-1; i++){
    if(photonPt[0] >= ptBins[i] && photonPt[0] < ptBins[i+1]){
      hVtxPtbinned[i] ->Fill(vtxN,weight*PUWeight);
      //hRho[i] -> Fill(rho,weight*PUWeight);
      break;
    }
    else if(photonPt[0]>=ptBins[numTrigger-1]){
      hVtxPtbinned[numTrigger-1]->Fill(vtxN,weight*PUWeight); 
      //hRho[numTrigger-1] -> Fill(rho,weight*PUWeight);
      break;
    }
  }

  
  hPhoton1Pt        -> Fill(photonPt[0],weight*PUWeight);
  hJet2Pt           -> Fill(jetPt2nd,weight*PUWeight);
  hJet1Pt           -> Fill(jetPt1stJet,weight*PUWeight);
  hMet              -> Fill(met_T1,weight*PUWeight);
  hMetPhiPhotonPhi  -> Fill(metPhi_T1,photonPhi[0],weight*PUWeight);
  hWeights          -> Fill(weight,PUWeight);
  hWeightsPhotonPt  -> Fill(weight,photonPt[0],PUWeight);

  hEcalRecHit      -> Fill(photonIsoEcal[0],weight*PUWeight);    
  hSigmaIetaIeta   -> Fill(photonSigmaIetaIeta[0],weight*PUWeight); 
  hHcalTowerSum    -> Fill(photonIsoHcal[0],weight*PUWeight); 
  hTrkHollowCone   -> Fill(photonIsoTrk[0],weight*PUWeight); 
  hHoverE          -> Fill(photonHadronicOverEM[0],weight*PUWeight); 
  pixelSeedVeto    -> Fill(photonHasPixelSeed[0],weight*PUWeight); 

  if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 1 || abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 2 || abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 3) hPhotonPtJetEtaFlavorLightQuarks_algo->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 4)  hPhotonPtJetEtaFlavorCharm_algo->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 5)  hPhotonPtJetEtaFlavorBottom_algo->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 21) hPhotonPtJetEtaFlavorGluon_algo->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_algo[corrJets.idx(idx1stJet)]) == 0)  hPhotonPtJetEtaFlavorNonDefined_algo->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);

  if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 1 || abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 2 || abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 3) hPhotonPtJetEtaFlavorLightQuarks_phys->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 4)  hPhotonPtJetEtaFlavorCharm_phys->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 5)  hPhotonPtJetEtaFlavorBottom_phys->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 21) hPhotonPtJetEtaFlavorGluon_phys->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);
  else if(abs(genJetID_phys[corrJets.idx(idx1stJet)]) == 0)  hPhotonPtJetEtaFlavorNonDefined_phys->Fill(photonPt[0],abs(jetEta[corrJets.idx(idx1stJet)]),weight*PUWeight);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  return 1;
  
}
