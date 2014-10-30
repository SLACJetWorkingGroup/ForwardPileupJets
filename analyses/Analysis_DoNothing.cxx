/**************************************************************************
 **
 **   File:         Analysis_DoNothing.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      M. Swiatlowski
 **
 **************************************************************************/

#define Analysis_DoNothing_cxx

#include "Analysis_DoNothing.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_DoNothing::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_DoNothing: DEBUG In WorkerBegin()" << endl;

  AnalysisBase::WorkerBegin();
  fDoLeptonSelection = false;
  ChainCfg()->Get("doLeptonSelection",  fDoLeptonSelection);


  if (Debug()) cout << "Analysis_DoNothing: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_DoNothing::ProcessEvent()
{

  if (Debug()) cout << "Analysis_DoNothing: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  cout << "Event Number " << Int("EventNumber") << endl;

 
  return true;
}


///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_DoNothing::WorkerTerminate()
{

  // Nothing

}
