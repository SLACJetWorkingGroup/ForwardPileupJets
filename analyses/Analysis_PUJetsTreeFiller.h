/**************************************************************************
 **
 **   File:         Analysis_PUJetsTreeFiller.h
 **
 **   Description:  Template for event-by-event tree filler
 **                 
 ** 
 **   Authors:      P. Nef
 **
 **   Created:      2/7/2014
 **
 **************************************************************************/

#ifndef Analysis_PUJetsTreeFiller_h
#define Analysis_PUJetsTreeFiller_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

 
using std::cout;
using std::endl;


class Analysis_PUJetsTreeFiller : public Analysis_JetMET_Base {

 public :
  
  Analysis_PUJetsTreeFiller(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_PUJetsTreeFiller() { }
  
  ClassDef(Analysis_PUJetsTreeFiller, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();


  TTree *fEventTree;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey);
  void FillEventVars(TTree *tree, const MomKey JetKey, Particle *myjet);
  void ResetBranches(TTree *tree);
  void MakeJetDisplay(Particle *myjet, const MomKey consttype);
  
  // per Event variables
  int   fTEventNumber;
  int   fTRunNumber;
  float fTEventWeight;
  float fTDefaultWeight;
  float fTMu;
  int   fTNPVtruth;
  int   fTNPV;
  float fTZvtx;

  // jets ----------------
  int   fTJetIndex;
  int   fTNJets;
  float fTJPt;
  float fTJEta;
  float fTJPhi;
  float fTJM;
  float fTJWidth;
  float fTdelRsqr;
  float fTClusSumPt;
  float fTdelRStdDev;
  float fTdelRSkewness;
  float fTdelRKurtosis;
  float fTdelR_01;
  float fTdelR_12;
  float fTdelR_23;
  float fTdelR_34;
  int  fTisHardScatter;
  int  fTisPileup;
  int  fTisQCDPileup1;
  int  fTisQCDPileup2;
  int  fTisQCDPileup3;
  int  fTisQCDPileup4;
  int  fTNumTowers;
  float fTTiming;
  float fTClusTimeAvg12;
  float fTClusTimeDif21;
  float fTClusTimeDif31;
  float fTClusTimeDif41;
  float fTClusTime1;
  float fTClusTime2;
  float fTClusTimeAdj;

  static const int MaxNCluster = 100;
  int fTClusIndex[MaxNCluster];
  float fTClusPt[MaxNCluster];
  float fTClusEta[MaxNCluster];
  float fTClusPhi[MaxNCluster];
  float fTcentlam[MaxNCluster];
  float fTEdens[MaxNCluster];
  float fTcellmaxfrac[MaxNCluster];
  float fTlong[MaxNCluster];
  float fTlat[MaxNCluster];
  float fTsecondLam[MaxNCluster];
  float fTsecondR[MaxNCluster];
  float fTiso[MaxNCluster];
  float fTsig[MaxNCluster];
  float fTEpos[MaxNCluster];
  float fTdelTheta[MaxNCluster];
  float fTdelPhi[MaxNCluster];
  float fTcentmag[MaxNCluster];
  float fTClusTime[MaxNCluster];
  int   fTNClus;
  int   fTIsLeadingClus[MaxNCluster];
  int   fTIs2LeadingClus[MaxNCluster];
  int   fTIs3LeadingClus[MaxNCluster];
  int   fTIs4LeadingClus[MaxNCluster];


    int fTNCaloTowers     ;
    float fTCaloTowersSumPt;
    float fTCaloTowersWidth;
    float fTCaloTowersWidthReCalc;
    float fTTopoTowersWidth;
    float fTTopoTowersWidthReCalc;

};

#endif

