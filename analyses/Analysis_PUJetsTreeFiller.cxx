/**************************************************************************
 **
 **   File:         Analysis_PUJetsTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PUJetsTreeFiller_cxx

#include "Analysis_PUJetsTreeFiller.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"
#include <TMath.h>

///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_PUJetsTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PUJetsTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();


  // trees -------------------
  fEventTree = new TTree("EventTree", "Tree with event-by-event variables");
  AddBranches(fEventTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PUJetsTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PUJetsTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  // Event Selection goes here... ---------------------------
  // only fill tree if event selection is good. this is set in Analysis_ForwardPileupJets.cxx
  if(Bool("EventSelection")==false) return true;  
  

  // Fill Tree-----------------------------------------------
  FillTree("AntiKt4LCTopoGood");


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_PUJetsTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PUJetsTreeFiller::WorkerTerminate()
{
    fEventTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_PUJetsTreeFiller::FillTree(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillTree Start" << endl;

  for(int iJ=0; iJ<jets(JetKey); ++iJ){
    Particle *myjet = &(jet(iJ, JetKey)); 
    if(myjet->p.Pt()>20){
      ResetBranches(fEventTree);
      fTNJets++;
      fTJetIndex = iJ;
      FillEventVars(fEventTree, JetKey, myjet);
      fEventTree->Fill();
//      MakeJetDisplay(myjet, "calotowersGhost");
//      MakeJetDisplay(myjet, "topotowersGhost");
    }
  
  }  
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillTree End" << endl;
  return;
}

///============================================
/// Make Event Display
///============================================
void Analysis_PUJetsTreeFiller::MakeJetDisplay( Particle *myjet, const MomKey consttype){
    
    if(myjet->p.Pt() < 20 || myjet->p.Pt()>30) return;
//    if(fabs(myjet->p.Eta())<2.4)               return;
//    if(myjet->Objs("constituents")!=1)         return;

    if(! myjet->Exists(consttype) )    return;

    cout << "MakeJetDisplay new jet: eta"  << myjet->p.Eta() << endl;
    ostringstream name;
    name << consttype << Int("EventNumber");
    float phi = myjet->p.Phi()>0? (round(10*(myjet->p.Phi()-0.01))/10.+0.02) : (round(10*(myjet->p.Phi()+0.01))/10.-0.02);
    float eta = myjet->p.Eta()>0? (round(10*(myjet->p.Eta()-0.05))/10.+0.05) : (round(10*(myjet->p.Eta()+0.05))/10.-0.05);

    TH2D* tower = new TH2D(name.str().c_str(), "", 11, eta-0.5-0.025, eta+0.5+0.025, 10, phi-0.5-0.025, phi+0.5+0.025);
    for(int iC = 0; iC < myjet->Objs(consttype); ++iC){
      Particle *constituent = (Particle*) myjet->Obj(consttype, iC);
      cout << "tower " << iC << " eta " << constituent->p.Eta() << " phi " << constituent->p.Phi() << endl;
      int bin = tower->FindBin(constituent->p.Eta(), constituent->p.Phi());
      tower->SetBinContent(bin, constituent->p.Pt());
    }
    ostringstream title; title << "1-cluster jet: ";
    if      (myjet->Int("isHardScatter")) title << "hard scatter";
    else if (myjet->Int("isPileup"))      title << "pileup";
    tower->SetTitle(title.str().c_str());
    tower->SetXTitle("#eta");
    tower->SetYTitle("#phi");
    tower->Write();
    delete tower;


    return;
}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_PUJetsTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("EventWeight" ,              &fTEventWeight,            "EventWeight/F");
    tree->Branch("DefaultWeight" ,            &fTDefaultWeight,          "DefaultWeight/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/I");
    tree->Branch("NPV" ,                      &fTNPV,                    "NPV/I");
    tree->Branch("Zvtx",                      &fTZvtx,                   "Zvtx/F");
  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("JetIndex",                  &fTJetIndex,               "JetIndex/I");
    tree->Branch("NJets",                     &fTNJets,                  "NJets/I");
    tree->Branch("Jpt",                       &fTJPt,                    "Jpt/F");
    tree->Branch("Jeta",                      &fTJEta,                   "Jeta/F");
    tree->Branch("Jphi",                      &fTJPhi,                   "Jphi/F");
    tree->Branch("Jm",                        &fTJM,                     "Jm/F");
    tree->Branch("Jwidth",                    &fTJWidth,                 "Jwidth/F");
    tree->Branch("delRsqr",                   &fTdelRsqr,                "delRsqr/F");
    tree->Branch("ClusSumPt",                 &fTClusSumPt,              "ClusSumPt/F");
    tree->Branch("delRStdDev",                &fTdelRStdDev,             "delRStdDev/F");
    tree->Branch("delRSkewness",              &fTdelRSkewness,           "delRSkewness/F");
    tree->Branch("delRKurtosis",              &fTdelRKurtosis,           "delRKurtosis/F");
    tree->Branch("delR_01",                   &fTdelR_01,                "delR_01/F");
    tree->Branch("delR_12",                   &fTdelR_12,                "delR_12/F");
    tree->Branch("delR_23",                   &fTdelR_23,                "delR_23/F");
    tree->Branch("delR_34",                   &fTdelR_34,                "delR_34/F");
    tree->Branch("isHardScatter",             &fTisHardScatter,          "isHardScatter/I");
    tree->Branch("isPileup",                  &fTisPileup,               "isPileup/I");
    tree->Branch("isQCDPileup1",              &fTisQCDPileup1,           "isQCDPileup1/I");
    tree->Branch("isQCDPileup2",              &fTisQCDPileup2,           "isQCDPileup2/I");
    tree->Branch("isQCDPileup3",              &fTisQCDPileup3,           "isQCDPileup3/I");
    tree->Branch("isQCDPileup4",              &fTisQCDPileup4,           "isQCDPileup4/I");
    tree->Branch("Timing",                    &fTTiming,                 "Timing/F");
    tree->Branch("NumTowers",                 &fTNumTowers,              "NumTowers/I");
    tree->Branch("NClus",                     &fTNClus,                  "NClus/I");
    tree->Branch("ClusTimeAvg12",             &fTClusTimeAvg12,          "fTClusTimeAvg12/F");
    tree->Branch("ClusTimeDif21",             &fTClusTimeDif21,          "fTClusTimeDif21/F");
    tree->Branch("ClusTimeDif31",             &fTClusTimeDif31,          "fTClusTimeDif31/F");
    tree->Branch("ClusTimeDif41",             &fTClusTimeDif41,          "fTClusTimeDif41/F");
    tree->Branch("ClusTime1",                 &fTClusTime1,              "ClusTime1/F");
    tree->Branch("ClusTime2",                 &fTClusTime2,              "ClusTime2/F");
    tree->Branch("ClusTimeAdj",               &fTClusTimeAdj,            "ClusTimeAdj/F");
    tree->Branch("ClusPt",                    &fTClusPt,                 "ClusPt[NClus]/F");
    tree->Branch("ClusEta",                   &fTClusEta,                "ClusEta[NClus]/F");
    tree->Branch("ClusPhi",                   &fTClusPhi,                "ClusPhi[NClus]/F");
    tree->Branch("centlam",                   &fTcentlam,                "centlam[NClus]/F");
    tree->Branch("Edens",                     &fTEdens,                  "Edens[NClus]/F");
    tree->Branch("cellmaxfrac",               &fTcellmaxfrac,            "cellmaxfrac[NClus]/F");
    tree->Branch("long",                      &fTlong,                   "long[NClus]/F");
    tree->Branch("lat",                       &fTlat,                    "lat[NClus]/F");
    tree->Branch("secondLam",                 &fTsecondLam,              "secondLam[NClus]/F");
    tree->Branch("secondR",                   &fTsecondR,                "secondR[NClus]/F");
    tree->Branch("iso",                       &fTiso,                    "iso[NClus]/F");
    tree->Branch("sig",                       &fTsig,                    "sig[NClus]/F");
    tree->Branch("Epos",                      &fTEpos,                   "Epos[NClus]/F");
    tree->Branch("delTheta",                  &fTdelTheta,               "delTheta[NClus]/F");
    tree->Branch("delPhi",                    &fTdelPhi,                 "delPhi[NClus]/F");
    tree->Branch("centmag",                   &fTcentmag,                "centmag[NClus]/F");
    tree->Branch("ClusIndex",                 &fTClusIndex,              "ClusIndex[NClus]/I");
    tree->Branch("ClusTime",                  &fTClusTime,               "ClusTime[NClus]/F");
    tree->Branch("ClusWidth",                 &fTClusWidth,              "ClusWidth[NClus]/F");
    tree->Branch("IsLeadingClus",             &fTIsLeadingClus,          "isLeadingClus[NClus]/I");
    tree->Branch("Is2LeadingClus",            &fTIs2LeadingClus,         "is2LeadingClus[NClus]/I");
    tree->Branch("Is3LeadingClus",            &fTIs3LeadingClus,         "is3LeadingClus[NClus]/I");
    tree->Branch("Is4LeadingClus",            &fTIs4LeadingClus,         "is4LeadingClus[NClus]/I");
    tree->Branch("NCaloTowers",               &fTNCaloTowers,            "NCaloTowers/I");
    tree->Branch("CaloTowersSumPt",           &fTCaloTowersSumPt,        "CaloTowersSumPt/F");
    tree->Branch("CaloTowersWidth",           &fTCaloTowersWidth,        "CaloTowersWidth/F");
    tree->Branch("CaloTowersWidthReCalc",     &fTCaloTowersWidthReCalc,  "CaloTowersWidthReCalc/F");
    tree->Branch("TopoTowersWidth",           &fTTopoTowersWidth,        "TopoTowersWidth/F");
    tree->Branch("TopoTowersWidthReCalc",     &fTTopoTowersWidthReCalc,  "TopoTowersWidthReCalc/F");

  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_PUJetsTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::ResetBranches Start" << endl;
  
    // Event Info
    fTEventNumber           = -999;
    fTRunNumber             = -999;
    fTEventWeight           = -999;
    fTDefaultWeight         = -999;
    fTNPVtruth              = -999;
    fTNPV                   = -999;
    fTMu                    = -999;
    fTZvtx                  = -999.99;
    
    //Jet Info
    fTJetIndex           = -999;
    fTJPt              = -999.99;
    fTJEta             = -999.99;
    fTJPhi             = -999.99;
    fTJM               = -999.99;
    fTJWidth           = -999.99;
    fTdelRsqr          = -999.99;
    fTClusSumPt        = -999.99;
    fTdelRStdDev       = -999.99;
    fTdelRSkewness     = -999.99;
    fTdelRKurtosis     = -999.99;
    fTdelR_01          = -999.99;
    fTdelR_12          = -999.99;
    fTdelR_23          = -999.99;
    fTdelR_34          = -999.99;
    fTisHardScatter    = -999;
    fTisPileup         = -999;
    fTisQCDPileup1     = -999;
    fTisQCDPileup2     = -999;
    fTisQCDPileup3     = -999;
    fTisQCDPileup4     = -999;
    fTNumTowers        = -999;
    fTTiming           = -999.99;
    fTNClus            = 0;
    fTClusTimeAvg12    = -999.99;
    fTClusTimeDif21    = -999.99;
    fTClusTimeDif31    = -999.99;
    fTClusTimeDif41    = -999.99;
    fTClusTime1        = -999.99;
    fTClusTime2        = -999.99;
    fTClusTimeAdj      = -999.99;
    for(int iC=0; iC<MaxNCluster; ++iC){
	fTClusPt[iC]     = -999.99;
        fTClusEta[iC]    = -999.99;
        fTClusPhi[iC]    = -999.99;
        fTcentlam[iC]    = -999.99;
        fTEdens[iC]      = -999.99;
        fTcellmaxfrac[iC]= -999.99;
	fTlong[iC]       = -999.99;
 	fTlat[iC]        = -999.99;
        fTsecondLam[iC]  = -999.99;
        fTsecondR[iC]    = -999.99;
	fTiso[iC]        = -999.99;
 	fTsig[iC]        = -999.99;
        fTEpos[iC]       = -999.99;
 	fTdelTheta[iC]   = -999.99;
        fTdelPhi[iC]     = -999.99;
        fTcentmag[iC]    = -999.99;
        fTClusTime[iC]       = -999.99;
        fTClusWidth[iC]       = -999.99;
        fTClusIndex[iC]  = -999;
        fTIsLeadingClus[iC] = -999;
        fTIs2LeadingClus[iC] = -999;
        fTIs3LeadingClus[iC] = -999;
        fTIs4LeadingClus[iC] = -999;
    }
    fTNCaloTowers     = 0;
    fTCaloTowersSumPt = -999;
    fTCaloTowersWidth = -999;
    fTCaloTowersWidthReCalc = -999;
    fTTopoTowersWidth = -999;
    fTTopoTowersWidthReCalc = -999;

  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_PUJetsTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey, Particle *myjet){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillEventVars Begin" << endl;
    
    // Event Info -----------------------------------
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTNPV                         = Exists("NPV")?      Int("NPV"):-1;
    fTEventWeight                 = Float("EventWeight");
    fTDefaultWeight               = DefaultWeight();
    fTMu                          = Float("averageIntPerXing");                  
    fTZvtx                        = vtx(0).x.z();

    // Jet Info ----------------------
    fTJPt         = myjet->p.Pt();
    fTJEta        = myjet->p.Eta();
    fTJPhi        = myjet->p.Phi();
    fTJM          = myjet->p.M();
    fTJWidth      = myjet->Float("WIDTH");
    fTdelRsqr     = myjet->Float("delRsqr");
    fTClusSumPt   = myjet->Float("ClusSumPt");
    fTdelRStdDev  = myjet->Float("delRStdDev");
    fTdelRSkewness = myjet->Float("delRSkewness");
    fTdelRKurtosis = myjet->Float("delRKurtosis"); 
    fTdelR_01     = myjet->Float("delR_01");
    fTdelR_12     = myjet->Float("delR_12");
    fTdelR_23     = myjet->Float("delR_23");
    fTdelR_34     = myjet->Float("delR_34");
    fTisHardScatter = myjet->Int("isHardScatter");
    fTisPileup      = myjet->Int("isPileup");
    if(isMC()){
      fTisQCDPileup1   = myjet->Int("isQCDPileup1");
      fTisQCDPileup2   = myjet->Int("isQCDPileup2");
      fTisQCDPileup3   = myjet->Int("isQCDPileup3");
      fTisQCDPileup4   = myjet->Int("isQCDPileup4");
    }
    fTNumTowers     = myjet->Float("NumTowers");   
    fTTiming     = myjet->Float("Timing");

    int LeadClusIndex = -999;
    int Lead2ClusIndex = -999;
    int Lead3ClusIndex = -999;
    int Lead4ClusIndex = -999;
    float LeadClusPt = 0;
    float Lead2ClusPt = 0;
    float Lead3ClusPt = 0;
    float Lead4ClusPt = 0;
    for(int iC=0; iC<myjet->Objs("constituents"); ++iC){
	Particle *con = (Particle*) myjet->Obj("constituents",iC);
	if(fTNClus == MaxNCluster) continue;
	fTClusIndex[iC] = iC;
        fTClusPt[iC] = con->p.Pt();
        fTClusEta[iC] = con->p.Eta();
        fTClusPhi[iC] = con->p.Phi();
        fTClusTime[iC] = con->Float("time");
        //if(fTClusPt[iC]>10){
          fTcentlam[iC] = con->Float("centerlambda");
          fTEdens[iC] = con->Float("firstEdens");
	  fTcellmaxfrac[iC] = con->Float("cellmaxfrac");
  	  fTlong[iC] = con->Float("longitudinal");
          fTlat[iC] = con->Float("lateral");
 	  fTsecondLam[iC] = con->Float("secondlambda");
	  fTsecondR[iC] = con->Float("secondR");
	  //fTiso[iC] = con->Float("isolation");
	  //fTsig[iC] = con->Float("significance");
	  //fTEpos[iC] = con->Float("eng_pos");
	  fTdelTheta[iC] = con->Float("deltaTheta");
	  fTdelPhi[iC] = con->Float("deltaPhi");
	  fTcentmag[iC] = con->Float("centermag");
	//}
	  fTClusWidth[iC] = 2.35*TMath::ATan(TMath::Sqrt(fTsecondR[iC])/fTcentmag[iC])*TMath::CosH(fTClusEta[iC]);
	  fTIsLeadingClus[iC] = 0;

	if(fTClusPt[iC]>LeadClusPt){
          LeadClusPt = fTClusPt[iC];
          LeadClusIndex = iC;
        }
        fTNClus++;
    }

    // Tower Information
    fTNCaloTowers     = myjet->Int("NCaloTowers");
    fTCaloTowersSumPt = myjet->Float("CaloTowersSumPt");
    fTCaloTowersWidth = myjet->Float("CaloTowersWidth");
    fTCaloTowersWidthReCalc = myjet->Float("CaloTowersWidthReCalc");
    fTTopoTowersWidth = myjet->Float("TopoTowersWidth");
    fTTopoTowersWidthReCalc = myjet->Float("TopoTowersWidthReCalc");

      for(int iC=0; iC<myjet->Objs("constituents"); ++iC){
        Particle *con = (Particle*) myjet->Obj("constituents",iC);
        fTClusPt[iC] = con->p.Pt();
        fTIs2LeadingClus[iC] = 0;
        if(iC==LeadClusIndex){continue;}
        if(fTClusPt[iC]<LeadClusPt && fTClusPt[iC]>Lead2ClusPt){
          Lead2ClusPt = fTClusPt[iC];
          Lead2ClusIndex = iC;
        }
      }
      for(int iC=0; iC<myjet->Objs("constituents"); ++iC){
        Particle *con = (Particle*) myjet->Obj("constituents",iC);
        fTClusPt[iC] = con->p.Pt();
        fTIs3LeadingClus[iC] = 0;
        if(iC==LeadClusIndex || iC==Lead2ClusIndex){continue;}
        if(fTClusPt[iC]<Lead2ClusPt && fTClusPt[iC]>Lead3ClusPt){
          Lead3ClusPt = fTClusPt[iC];
          Lead3ClusIndex = iC;
        }
      }
      for(int iC=0; iC<myjet->Objs("constituents"); ++iC){
        Particle *con = (Particle*) myjet->Obj("constituents",iC);
        fTClusPt[iC] = con->p.Pt();
        fTIs4LeadingClus[iC] = 0;
        if(iC==LeadClusIndex || iC==Lead2ClusIndex || iC==Lead3ClusIndex){continue;}
        if(fTClusPt[iC]<Lead3ClusPt && fTClusPt[iC]>Lead4ClusPt){
          Lead4ClusPt = fTClusPt[iC];
          Lead4ClusIndex = iC;
        }
      }
    
   fTIsLeadingClus[LeadClusIndex] = 1;
   fTIs2LeadingClus[Lead2ClusIndex] = 1;
   fTIs3LeadingClus[Lead3ClusIndex] = 1;
   fTIs4LeadingClus[Lead4ClusIndex] = 1;

  if(myjet->Objs("constituents")<2){
    fTIs2LeadingClus[Lead2ClusIndex] = 0;    
    fTIs3LeadingClus[Lead3ClusIndex] = 0;
    fTIs4LeadingClus[Lead4ClusIndex] = 0;
  }
  if(myjet->Objs("constituents")<3){
    fTIs3LeadingClus[Lead3ClusIndex] = 0;
    fTIs4LeadingClus[Lead4ClusIndex] = 0;
  }
  if(myjet->Objs("constituents")<4){
    fTIs4LeadingClus[Lead4ClusIndex] = 0;
  }
  fTClusTimeAvg12 = (fTClusTime[LeadClusIndex]+fTClusTime[Lead2ClusIndex])/2;
  fTClusTimeDif21 = fTClusTime[Lead2ClusIndex]-fTClusTime[LeadClusIndex];
  fTClusTimeDif31 = fTClusTime[Lead3ClusIndex]-fTClusTime[LeadClusIndex];
  fTClusTimeDif41 = fTClusTime[Lead4ClusIndex]-fTClusTime[LeadClusIndex];
  fTClusTime1 = fTClusTime[LeadClusIndex];
  fTClusTime2 = fTClusTime[Lead2ClusIndex];

  TFile f("~/nfs/ProofAna20140926/ForwardPileupJets/analyses/TimeCorrection.root");
  TH2D *TimeMap = (TH2D*)f.Get("TimeMap");
  fTClusTimeAdj = fTClusTime[LeadClusIndex]-TimeMap->GetBinContent(TimeMap->FindBin(fTZvtx,fTJEta));

  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillEventVars End" << endl;
  return;
}


