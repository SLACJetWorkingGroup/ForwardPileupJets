/**************************************************************************
 **
 **   File:         Analysis_ForwardPileupJets.h
 **
 **   Description:  Computation and plotting of ForwardPileupJets studies.
 **                 
 ** 
 **   Authors:      M. Swiatlowski
 **
 **   Created:      2013-01-22
 **
 **************************************************************************/

#ifndef Analysis_ForwardPileupJets_h
#define Analysis_ForwardPileupJets_h

#include "Analysis_ForwardPileupJets_Base.h"
 
using std::cout;
using std::endl;



class Analysis_ForwardPileupJets : public Analysis_ForwardPileupJets_Base {

 public :
  
  Analysis_ForwardPileupJets(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_ForwardPileupJets() { }
  
  ClassDef(Analysis_ForwardPileupJets, 0);
  
  
  
  Bool_t  fDetail;
  Bool_t  fDoLeptonSelection;
  Bool_t  fDoQCDSelection;
  Bool_t  fDoTowers;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();
  
  
  void   CalculateTowerJetMoments(MomKey JetKey);
  
  bool   EventSelection();
  bool   MakeJetMassCut(MomKey JetKey);
  bool   MakeJetPtCut(MomKey JetKey);
  void   MakeJetCollections(const MomKey JetKey);

  void   CalculateVolatilities(const MomKey JetKey, const MomKey QKey, const int NJetsCompute = 2);

  // helper functions to facilitate similar types of plots
  void MakeJetPlots(const MomKey JetKey1, const MomKey JetKey2, const MomKey JetKey3);

  float GetConstitSumPt(Particle *myjet, const MomKey constType="constituents");
  float GetConstitPtWeightedMeanOverDr(Particle *myjet, const MomKey constType="constituents");
  float GetConstitPtWeightedStdDevOverDr(Particle *myjet, const MomKey constType="constituents");
  float GetConstitPtWeightedSkewnessOverDr(Particle *myjet, const MomKey constType="constituents");
  float GetConstitPtWeightedKurtosisOverDr(Particle *myjet, const MomKey constType="constituents");

  private :			  

};

#endif

