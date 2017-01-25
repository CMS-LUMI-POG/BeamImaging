/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "gauss2dProjectionV1_Ext.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
#include "TMatrixDSym.h"
#include "TMatrixD.h"

ClassImp(gauss2dProjectionV1_Ext) 

 gauss2dProjectionV1_Ext::gauss2dProjectionV1_Ext(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _y,
                        RooAbsReal& _x01,
                        RooAbsReal& _y01,
                        RooAbsReal& _x_Width1,
                        RooAbsReal& _y_Width1,
                        RooAbsReal& _rho1,
                        RooAbsReal& _x02,
                        RooAbsReal& _y02,
                        RooAbsReal& _x_Width2,
                        RooAbsReal& _y_Width2,
                        RooAbsReal& _rho2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   y("y","y",this,_y),
   x01("x01","x01",this,_x01),
   y01("y01","y01",this,_y01),
   x_Width1("x_Width1","x_Width1",this,_x_Width1),
   y_Width1("y_Width1","y_Width1",this,_y_Width1),
   rho1("rho1","rho1",this,_rho1),
   x02("x02","x02",this,_x02),
   y02("y02","y02",this,_y02),
   x_Width2("x_Width2","x_Width2",this,_x_Width2),
   y_Width2("y_Width2","y_Width2",this,_y_Width2),
   rho2("rho2","rho2",this,_rho2)
 { 
 } 


 gauss2dProjectionV1_Ext::gauss2dProjectionV1_Ext(const gauss2dProjectionV1_Ext& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   y("y",this,other.y),
   x01("x01",this,other.x01),
   y01("y01",this,other.y01),
   x_Width1("x_Width1",this,other.x_Width1),
   y_Width1("y_Width1",this,other.y_Width1),
   rho1("rho1",this,other.rho1),
   x02("x02",this,other.x02),
   y02("y02",this,other.y02),
   x_Width2("x_Width2",this,other.x_Width2),
   y_Width2("y_Width2",this,other.y_Width2),
   rho2("rho2",this,other.rho2)
 { 
 } 



 Double_t gauss2dProjectionV1_Ext::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 

TMatrixDSym beam1(2);
    beam1(0,0) = TMath::Power(x_Width1,2.0); 
    beam1(1,1) = TMath::Power(y_Width1,2.0); 
    beam1(1,0) = rho1*x_Width1*y_Width1;
    beam1(0,1) = rho1*x_Width1*y_Width1;
    TMatrixDSym beam2MargInv(2);
    beam2MargInv(0,0) = 0.;
    beam2MargInv(1,1) =  1./TMath::Power(y_Width2,2.0);;
    beam2MargInv(0,1) = 0.;
    beam2MargInv(1,0) = 0.;
    TMatrixDSym vtxRes(2);
    vtxRes(0,0) = 0.49;
    vtxRes(1,1) = 0.49;
    vtxRes(0,1) = 0.;
    vtxRes(1,0) = 0.;
    TMatrixDSym sigmaFinalInv(2);
    sigmaFinalInv = ((beam1.Invert()+beam2MargInv).Invert()+vtxRes).Invert();
    fitFunc->SetParameter(0,sigmaFinalInv(0,0));
    fitFunc->SetParameter(1,sigmaFinalInv(1,1));
    fitFunc->SetParameter(2,sigmaFinalInv(1,0));
    fitFunc->SetParameter(3,sigmaFinalInv.Invert().Determinant());

   // return n2*(twodGaussProj->Integral(-15,15,y,y+0.00000001))* n1*twodGauss->Eval(x,y); 

    TMatrixDSym sigmaFinal(2);
    TMatrixD b10(2,1);
    TMatrixD b20(2,1);
    b10(0,0) = x01;
    b10(1,0) = y01;
    sigmaFinal = (sigmaFinalInv.Invert()-vtxRes);
    b20(0,0) = x02;
    b20(1,0) = y02;
    TMatrixD c0 = (sigmaFinal*beam1.Invert())*b10 + (sigmaFinal*beam2MargInv)*b20;

    //std::cout<<"C0 "<<c0(0,0)<<"  "<<c0(1,0)<<std::endl;
    //std::cout<<"B10 before"<<x01<<"  "<<y01<<std::endl; 
    //std::cout<<"B10 "<<b10(0,0)<<"  "<<b10(1,0)<<std::endl;
    
   //return twodGauss1->Eval(x,y);

    return fitFunc->Eval((x-x01),(y-y01));

    //return fitFunc->Eval((x-c0(0,0)),(y-c0(1,0)));

 } 


