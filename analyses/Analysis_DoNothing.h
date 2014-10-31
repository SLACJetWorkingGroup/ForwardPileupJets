/**************************************************************************
 **
 **   File:         Analysis_DoNothing.h
 **
 **   Description:  Computation and plotting of DoNothing studies.
 **                 
 ** 
 **   Authors:      M. Swiatlowski
 **
 **   Created:      2013-01-22
 **
 **************************************************************************/

#ifndef Analysis_DoNothing_h
#define Analysis_DoNothing_h

#include "AnalysisBase.h"
 
using std::cout;
using std::endl;



class Analysis_DoNothing : public AnalysisBase {

 public :
  
  Analysis_DoNothing(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_DoNothing() { }
  
  ClassDef(Analysis_DoNothing, 0);
  
  
  
  Bool_t  fDetail;
  Bool_t  fDoLeptonSelection;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();
  
  
  
  bool   EventSelection();
  bool   MakeJetMassCut(MomKey JetKey);
  bool   MakeJetPtCut(MomKey JetKey);
  void   MakeJetCollections(const MomKey JetKey);

  void   CalculateVolatilities(const MomKey JetKey, const MomKey QKey, const int NJetsCompute = 2);

  // helper functions to facilitate similar types of plots
  void MakeJetPlots(const MomKey JetKey1, const MomKey JetKey2, const MomKey JetKey3);

  float GetConstitSumPt(Particle *myjet);
  float GetConstitPtWeightedMeanOverDr(Particle *myjet);
  float GetConstitPtWeightedStdDevOverDr(Particle *myjet);
  float GetConstitPtWeightedSkewnessOverDr(Particle *myjet);
  float GetConstitPtWeightedKurtosisOverDr(Particle *myjet);

  private :			  

};

#endif

