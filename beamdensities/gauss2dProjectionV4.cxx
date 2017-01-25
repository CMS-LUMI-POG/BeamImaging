/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "gauss2dProjectionV4.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(gauss2dProjectionV4) 

 gauss2dProjectionV4::gauss2dProjectionV4(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _y,
                        RooAbsReal& _n1,
                        RooAbsReal& _x_Width1,
                        RooAbsReal& _y_Width1,
                        RooAbsReal& _rho1,
                        RooAbsReal& _n2,
                        RooAbsReal& _x_Width2,
                        RooAbsReal& _y_Width2,
                        RooAbsReal& _rho2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   y("y","y",this,_y),
   n1("n1","n1",this,_n1),
   x_Width1("x_Width1","x_Width1",this,_x_Width1),
   y_Width1("y_Width1","y_Width1",this,_y_Width1),
   rho1("rho1","rho1",this,_rho1),
   n2("n2","n2",this,_n2),
   x_Width2("x_Width2","x_Width2",this,_x_Width2),
   y_Width2("y_Width2","y_Width2",this,_y_Width2),
   rho2("rho2","rho2",this,_rho2)
 { 
 } 


 gauss2dProjectionV4::gauss2dProjectionV4(const gauss2dProjectionV4& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   y("y",this,other.y),
   n1("n1",this,other.n1),
   x_Width1("x_Width1",this,other.x_Width1),
   y_Width1("y_Width1",this,other.y_Width1),
   rho1("rho1",this,other.rho1),
   n2("n2",this,other.n2),
   x_Width2("x_Width2",this,other.x_Width2),
   y_Width2("y_Width2",this,other.y_Width2),
   rho2("rho2",this,other.rho2)
 { 
 } 



 Double_t gauss2dProjectionV4::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
 //TF2 *twodGauss = new TF2("twodGauss","exp(-0.5*((x[0])/[0])^2-0.5*TMath::Power(((x[1])/[1]),2)-2*[2]*(x[0])*(x[1])/([0]*[1]))",-15,15,-15,15);
   twodGauss->SetParameter(0,x_Width1);
   twodGauss->SetParameter(1,y_Width1);
   twodGauss->SetParameter(2,rho1);

   //TF2 *twodGaussProj = new TF2("twodGaussProj","exp(-0.5*((x[0])/[0])^2-0.5*TMath::Power(((x[1])/[1]),2)-2*[2]*(x[0])*(x[1])/([0]*[1]))",-15,15,-15,15);
   twodGaussProj->SetParameter(0,x_Width2);
   twodGaussProj->SetParameter(1,y_Width2);
   twodGaussProj->SetParameter(2,rho2);
   
   double convVal;

   for(double dx = -15; dx<15; dx+=epsilon){
     for(double dy = -15; dy<15; dy+=epsilon){
       convVal+=n2*(projectTF2(x-dx,y-dy)*twodGauss->Eval(x-dx,y-dy))*n1*vtxResolution->Eval(dx,dy)*epsilon*epsilon;
     }
   }
   std::cout<<convVal<<std::endl;
   return  convVal ; 
 } 


