/**************************************************************************
 **
 **   File:         Analysis_TopCommonSelectionForward.h
 **
 **   Description:  Analysis class for common selection of top-like events
 **                 
 ** 
 **   Authors:      M. Swiatlowski, B. Nachman
 **
 **   Created:      2012-01-30
 **
 **************************************************************************/

#ifndef Analysis_TopCommonSelectionForward_h
#define Analysis_TopCommonSelectionForward_h

#include "Analysis_JetMET_Base.h"
 
using std::cout;
using std::endl;

namespace Root{
  class TTileTripReader;
}


class Analysis_TopCommonSelectionForward : public Analysis_JetMET_Base {

 public :
  
  Analysis_TopCommonSelectionForward(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }

  virtual ~Analysis_TopCommonSelectionForward() { }
  
  ClassDef(Analysis_TopCommonSelectionForward, 0);
  
  Bool_t  fDetail;

  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();
  void    Book(const TString& prefix);
  void    BookBasic();
  void    BookDetail();
  

 private :			  
  bool doSMWZ;
  bool doCOMMON;
  int SystType;
  int doNOJVFcut;

  MomKey BKey;
  MomKey MuonSFKey;
  MomKey MuonTrigSFKey;

  Root::TTileTripReader* m_treader;

};

#endif

