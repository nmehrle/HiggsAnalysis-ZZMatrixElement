#ifndef _TEVTPROB_HH_
#define _TEVTPROB_HH_
//-----------------------------------------------------------------------------
// Description: Class TEvtProb: EvtProb base class
// ------------
//
//      Event Probability Density Calculation
//
// Feb 21 2011
// Sergo Jindariani
// Yanyan Gao
//-----------------------------------------------------------------------------
#include <sstream>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "assert.h"
#include "TROOT.h"
// ME related
#include "TMCFM.hh"
#include "TVar.hh"
#include "TUtil.hh"
//#include "Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h"
#include "ZZMatrixElement/MELA/interface/HiggsCSandWidth_MELA.h"


//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject {
  
public:
  //--------------------
  // Variables
  //--------------------
  TVar::Process _process;
  TVar::MatrixElement _matrixElement;
  TVar::Production _production;
  HiggsCSandWidth_MELA *myCSW_;
  //HiggsCSandWidth *myCSW_;
  double _hmass;
  double _hwidth;
  double EBEAM;
  
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb() {};
  TEvtProb(const char* path,double ebeam);
 // TEvtProb(const string* path,double ebeam);
  ~TEvtProb();
  
  //----------------------
  // Function
  //----------------------
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void SetProduction(TVar::Production tmp){ _production = tmp; }
  double XsecCalc(TVar::Process proc,
		  TVar::Production production,
		  const hzz4l_event_type &hzz4l_event,
		  TVar::VerbosityLevel verbosity, double couplingvals[2], double selfDHvvcoupl[30][2],double selfDZqqcoupl[2][2],
     double selfDZvvcoupl[2][2],
     double selfDGqqcoupl[2][2],
     double selfDGggcoupl[5][2],
     double selfDGvvcoupl[10][2]);
  double XsecCalcXJJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[3],
		     TVar::VerbosityLevel verbosity);

  // this appears to be some kind of 
  // way of setting MCFM parameters through
  // an interface defined in TMCFM.hh
  void SetHiggsMass(double mass, float wHiggs=-1);
  ClassDef(TEvtProb,0);
};

#endif

