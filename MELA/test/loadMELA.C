//
// MELA root package loader - see testKD.C for instructions
//
{
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->Load("libZZMatrixElementMELA.so");
//  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/libmcfm.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h+");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/TVar.hh+");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/src/computeAngles.h+");
//  gROOT->LoadMacro("RooSpinZero_KD_ZH.cc+");
//  gROOT->LoadMacro("HZZ4L_RooSpinZeroPdf.cc+");
  //gSystem->AddIncludePath("-I$CMSSW_BASE/src/ ");  
}
