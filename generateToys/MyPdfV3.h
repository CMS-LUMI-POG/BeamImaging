/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef MYPDFV3
#define MYPDFV3

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TF2.h"
 
class MyPdfV3 : public RooAbsPdf {
public:
  MyPdfV3() {} ; 
  MyPdfV3(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _y,
	      RooAbsReal& _wNW1,
	      RooAbsReal& __rhoN1,
	      RooAbsReal& _sigxN1,
	      RooAbsReal& _sigyN1,
	      RooAbsReal& __rhoW1,
	      RooAbsReal& _sigxW1,
	      RooAbsReal& _sigyW1,
	      RooAbsReal& _wNW2,
	      RooAbsReal& _sigyN2,
	      RooAbsReal& _sigyW2);
  MyPdfV3(const MyPdfV3& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new MyPdfV3(*this,newname); }
  inline virtual ~MyPdfV3() { }

protected:

  RooRealProxy x ;
  RooRealProxy y ;
  RooRealProxy wNW1 ;
  RooRealProxy _rhoN1 ;
  RooRealProxy sigxN1 ;
  RooRealProxy sigyN1 ;
  RooRealProxy _rhoW1 ;
  RooRealProxy sigxW1 ;
  RooRealProxy sigyW1 ;
  RooRealProxy wNW2 ;
  RooRealProxy sigyN2 ;
  RooRealProxy sigyW2 ;
  
  TF2 *fitFuncN1N2 = new TF2("fitFuncN1N2","[4]*1./(2*3.14159*TMath::Sqrt([3]))*TMath::Exp(-0.5*(TMath::Power(x[0],2.0)*[0]+TMath::Power(x[1],2.0)*[1]+2*[2]*x[0]*x[1]))",-30,30,-30,30);
  TF2 *fitFuncN1W2 = new TF2("fitFuncN1W2","[4]*1./(2*3.14159*TMath::Sqrt([3]))*TMath::Exp(-0.5*(TMath::Power(x[0],2.0)*[0]+TMath::Power(x[1],2.0)*[1]+2*[2]*x[0]*x[1]))",-30,30,-30,30);
  TF2 *fitFuncW1N2 = new TF2("fitFuncW1N2","[4]*1./(2*3.14159*TMath::Sqrt([3]))*TMath::Exp(-0.5*(TMath::Power(x[0],2.0)*[0]+TMath::Power(x[1],2.0)*[1]+2*[2]*x[0]*x[1]))",-30,30,-30,30);
  TF2 *fitFuncW1W2 = new TF2("fitFuncW1W2","[4]*1./(2*3.14159*TMath::Sqrt([3]))*TMath::Exp(-0.5*(TMath::Power(x[0],2.0)*[0]+TMath::Power(x[1],2.0)*[1]+2*[2]*x[0]*x[1]))",-30,30,-30,30);

  Double_t evaluate() const ;

private:

  ClassDef(MyPdfV3,1) // Your description goes here...
};
 
#endif
