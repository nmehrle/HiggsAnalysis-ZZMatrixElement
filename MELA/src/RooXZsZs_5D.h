/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOXZSZS_5D
#define ROOXZSZS_5D

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooXZsZs_5D : public RooAbsPdf {
public:
    RooXZsZs_5D() {} ; 
    RooXZsZs_5D(const char *name, const char *title,
                RooAbsReal& _m1,
                RooAbsReal& _m2,
                RooAbsReal& _h1,
                RooAbsReal& _h2,
                RooAbsReal& _Phi,
                RooAbsReal& _a1Val,
                RooAbsReal& _phi1Val,
                RooAbsReal& _a2Val,
                RooAbsReal& _phi2Val,
                RooAbsReal& _a3Val,
                RooAbsReal& _phi3Val,
		RooAbsReal& _useGTerm,
		RooAbsReal& _g1Val,
		RooAbsReal& _g2Val,
		RooAbsReal& _g3Val,
		RooAbsReal& _g4Val,
    RooAbsReal& _g1_primeVal,
    RooAbsReal& _g2_primeVal,
    RooAbsReal& _g3_primeVal,
    RooAbsReal& _g4_primeVal,
    RooAbsReal& _g1_prime2Val,
    RooAbsReal& _g2_prime2Val,
    RooAbsReal& _g3_prime2Val,
    RooAbsReal& _g4_prime2Val,
    RooAbsReal& _g1_prime3Val,
    RooAbsReal& _g2_prime3Val,
    RooAbsReal& _g3_prime3Val,
    RooAbsReal& _g4_prime3Val,
    RooAbsReal& _g1_prime4Val,
    RooAbsReal& _g2_prime4Val,
    RooAbsReal& _g3_prime4Val,
    RooAbsReal& _g4_prime4Val,
                RooAbsReal& _mZ,
                RooAbsReal& _gamZ,
                RooAbsReal& _mX,
                RooAbsReal& _R1Val,
                RooAbsReal& _R2Val);
    RooXZsZs_5D(const RooXZsZs_5D& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const { return new RooXZsZs_5D(*this,newname); }
    inline virtual ~RooXZsZs_5D() { }
    
    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
    
protected:
    
    RooRealProxy m1 ;
    RooRealProxy m2 ;
    RooRealProxy h1 ;
    RooRealProxy h2 ;
    RooRealProxy Phi ;
    RooRealProxy a1Val ;
    RooRealProxy phi1Val ;
    RooRealProxy a2Val ;
    RooRealProxy phi2Val ;
    RooRealProxy a3Val ;
    RooRealProxy phi3Val ;
    RooRealProxy useGTerm ;
    RooRealProxy g1Val ;
    RooRealProxy g2Val ;
    RooRealProxy g3Val ;
    RooRealProxy g4Val ;
    RooRealProxy g1_primeVal;
    RooRealProxy g2_primeVal;
    RooRealProxy g3_primeVal;
    RooRealProxy g4_primeVal;
    RooRealProxy g1_prime2Val;
    RooRealProxy g2_prime2Val;
    RooRealProxy g3_prime2Val;
    RooRealProxy g4_prime2Val;
    RooRealProxy g1_prime3Val;
    RooRealProxy g2_prime3Val;
    RooRealProxy g3_prime3Val;
    RooRealProxy g4_prime3Val;
    RooRealProxy g1_prime4Val;
    RooRealProxy g2_prime4Val;
    RooRealProxy g3_prime4Val;
    RooRealProxy g4_prime4Val;
    RooRealProxy mZ ;
    RooRealProxy gamZ ;
    RooRealProxy mX ;
    RooRealProxy R1Val ;
    RooRealProxy R2Val ;
    
    Double_t evaluate() const ;
    
private:
    double Lambda_z1 ; 
    double Lambda_z2 ; 
    double Lambda_z3 ; 
    double Lambda_z4 ; 
    //ClassDef(RooXZsZs_5D,1) // Your description goes here...
};

#endif
