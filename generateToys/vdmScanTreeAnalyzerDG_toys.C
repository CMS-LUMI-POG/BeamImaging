#include "TH2F.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TMath.h"
#include "TString.h"
//#include "RooFit.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <sstream>
#include "TRint.h"
//#include "TROOT.h"
//using namespace std;
//#include "RooWorkspace.h"
#include "TTree.h"
//gSystem->Load("libRooFit");
//using namespace RooFit ;
//#include "TMatrixDSym.h"
#include "MyPdfV3.h"
#include "MyPdfV4.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooCFunction1Binding.h" 
#include "RooCFunction2Binding.h"
#include "RooTFnBinding.h"
#include "RooMultiVarGaussian.h"
#include "RooProdPdf.h"
#include "RooNumConvPdf.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
using namespace RooFit;


Double_t beamMultDG(Double_t *x,Double_t *par)
{
  Double_t arg = 0;
  Double_t pi = 3.1415926;

  double x01 = par[0];
  double y01 = par[1];
  double xwidthN1 = par[2];
  double ywidthN1 = par[3];
  double rhoN1 = par[4];
  double xwidthW1 = par[5];
  double ywidthW1 = par[6];
  double rhoW1 = par[7];
  double nw_weight1 = par[8];

  double x02 = par[9];
  double y02 = par[10];
  double xwidthN2 = par[11];
  double ywidthN2 = par[12];
  double rhoN2 = par[13];
  double xwidthW2 = par[14];
  double ywidthW2 = par[15];
  double rhoW2 = par[16];
  double nw_weight2 = par[17];
  
  double beamN1_ = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN1,2)))*TMath::Sqrt(TMath::Power(xwidthN1,2.)/*+0.49*/)*TMath::Sqrt(TMath::Power(ywidthN1,2.)/*+0.49*/))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN1,2))*(TMath::Power((x[0]-x01)/xwidthN1,2.0)+TMath::Power((x[1]-y01)/ywidthN1,2.0)-2*rhoN1*(x[0]-x01)*(x[1]-y01)/(xwidthN1*ywidthN1)));

  double beamW1_ = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW1,2)))*TMath::Sqrt(TMath::Power(xwidthW1,2.)/*+0.49*/)*TMath::Sqrt(TMath::Power(ywidthW1,2.)/*+0.49*/))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW1,2))*(TMath::Power((x[0]-x01)/xwidthW1,2.0)+TMath::Power((x[1]-y01)/ywidthW1,2.0)-2*rhoW1*(x[0]-x01)*(x[1]-y01)/(xwidthW1*ywidthW1)));

  double beamN2_ = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN2,2)))*TMath::Sqrt(TMath::Power(xwidthN2,2.)/*+0.49*/)*TMath::Sqrt(TMath::Power(ywidthN2,2.)/*+0.49*/))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN2,2))*(TMath::Power((x[0]-x02)/xwidthN2,2.0)+TMath::Power((x[1]-y02)/ywidthN2,2.0)-2*rhoN2*(x[0]-x02)*(x[1]-y02)/(xwidthN2*ywidthN2)));

  double beamW2_ = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW2,2)))*TMath::Sqrt(TMath::Power(xwidthW2,2.)/*+0.49*/)*TMath::Sqrt(TMath::Power(ywidthW2,2.)/*+0.49*/))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW2,2))*(TMath::Power((x[0]-x02)/xwidthW2,2.0)+TMath::Power((x[1]-y02)/ywidthW2,2.0)-2*rhoW2*(x[0]-x02)*(x[1]-y02)/(xwidthW2*ywidthW2)));

  double product2 = (nw_weight1*beamN1_ + (1.-nw_weight1)*beamW1_) * (nw_weight2*beamN2_ + (1.-nw_weight2)*beamW2_);

  return product2;
 
}



Double_t beamMultDGt(Double_t *x,Double_t *par)
{
  Double_t arg = 0;
  Double_t pi = 3.1415926;

  double x01 = par[0];
  double y01 = par[1];
  double xwidthN1 = par[2];
  double ywidthN1 = par[3];
  double rhoN1 = par[4];
  double xwidthW1 = par[5];
  double ywidthW1 = par[6];
  double rhoW1 = par[7];
  double nw_weight1 = par[8];

  double x02 = par[9];
  double y02 = par[10];
  double xwidthN2 = par[11];
  double ywidthN2 = par[12];
  double rhoN2 = par[13];
  double xwidthW2 = par[14];
  double ywidthW2 = par[15];
  double rhoW2 = par[16];
  double nw_weight2 = par[17];

TMatrixDSym sigN1(2);
  sigN1(0,0) = TMath::Power(xwidthN1,2.);
  sigN1(1,1) = TMath::Power(ywidthN1,2.);
  sigN1(1,0) = rhoN1*xwidthN1*ywidthN1;
  sigN1(0,1) = rhoN1*xwidthN1*ywidthN1;
  double sigN1_det = sigN1.Determinant();
  sigN1.Invert();
  //sigN1.Print(); 
  TMatrixDSym sigW1(2);
  sigW1(0,0) = TMath::Power(xwidthW1,2.);
  sigW1(1,1) = TMath::Power(ywidthW1,2.);
  sigW1(1,0) = rhoW1*xwidthW1*ywidthW1;
  sigW1(0,1) = rhoW1*xwidthW1*ywidthW1;
  double sigW1_det = sigW1.Determinant();
  sigW1.Invert();
  //sigW1.Print();
  TMatrixDSym sigN2(2);
  sigN2(0,0) = TMath::Power(xwidthN2,2.);
  sigN2(1,1) = TMath::Power(ywidthN2,2.);
  sigN2(1,0) = rhoN2*xwidthN2*ywidthN2;
  sigN2(0,1) = rhoN2*xwidthN2*ywidthN2;
  double sigN2_det = sigN2.Determinant();
  sigN2.Invert();
  //sigN2.Print();
  TMatrixDSym sigW2(2);
  sigW2(0,0) = TMath::Power(xwidthW2,2.);
  sigW2(1,1) = TMath::Power(ywidthW2,2.);
  sigW2(1,0) = rhoW2*xwidthW2*ywidthW2;
  sigW2(0,1) = rhoW2*xwidthW2*ywidthW2;
  double sigW2_det = sigW2.Determinant();
  sigW2.Invert();
  //sigW2.Print();

  double beamN1 = (nw_weight1)*1./(TMath::Sqrt(sigN1_det)*2*pi)*exp(-0.5*(TMath::Power((x[0]-x01),2.0)*sigN1(0,0)+TMath::Power(x[1]-y01,2.0)*sigN1(1,1)+2*sigN1(1,0)*(x[0]-x01)*(x[1]-y01)));
  double beamN2 = (nw_weight2)*1./(TMath::Sqrt(sigN2_det)*2*pi)*exp(-0.5*(TMath::Power((x[0]-x02),2.0)*sigN2(0,0)+TMath::Power(x[1]-y02,2.0)*sigN2(1,1)+2*sigN2(1,0)*(x[0]-x02)*(x[1]-y02)));

  double beamW1 = (1-nw_weight1)*1./(TMath::Sqrt(sigW1_det)*2*pi)*exp(-0.5*(TMath::Power((x[0]-x01),2.0)*sigW1(0,0)+TMath::Power(x[1]-y01,2.0)*sigW1(1,1)+2*sigW1(1,0)*(x[0]-x01)*(x[1]-y01)));

  double beamW2 = (1-nw_weight2)*1./(TMath::Sqrt(sigW2_det)*2*pi)*exp(-0.5*(TMath::Power((x[0]-x02),2.0)*sigW2(0,0)+TMath::Power(x[1]-y02,2.0)*sigW2(1,1)+2*sigW2(1,0)*(x[0]-x02)*(x[1]-y02)));
  
 
  


  /* TMatrixDSym sigN1(2);
  sigN1(0,0) = TMath::Power(xwidthN1,2.);
  sigN1(1,1) = TMath::Power(ywidthN1,2.);
  sigN1(1,0) = rhoN1*xwidthN1*ywidthN1;
  sigN1(0,1) = rhoN1*xwidthN1*ywidthN1;
  double sigN1_det = sigN1.Determinant();
  sigN1.Invert();
  TMatrixDSym sigW1(2);
  sigW1(0,0) = TMath::Power(xwidthW1,2.);
  sigW1(1,1) = TMath::Power(ywidthW1,2.);
  sigW1(1,0) = rhoW1*xwidthW1*ywidthW1;
  sigW1(0,1) = rhoW1*xwidthW1*ywidthW1;
  double sigW1_det = sigW1.Determinant();
  sigW1.Invert();
  TMatrixDSym sigN2(2);
  sigN2(0,0) = TMath::Power(xwidthN2,2.);
  sigN2(1,1) = TMath::Power(ywidthN2,2.);
  sigN2(1,0) = rhoN2*xwidthN2*ywidthN2;
  sigN2(0,1) = rhoN2*xwidthN2*ywidthN2;
  double sigN2_det = sigN2.Determinant();
  sigN2.Invert();
  TMatrixDSym sigW2(2);
  sigW2(0,0) = TMath::Power(xwidthW2,2.);
  sigW2(1,1) = TMath::Power(ywidthW2,2.);
  sigW2(1,0) = rhoW2*xwidthW2*ywidthW2;
  sigW2(0,1) = rhoW2*xwidthW2*ywidthW2;
  double sigW2_det = sigW2.Determinant();
  sigW2.Invert();

  double beamN1 = (nw_weight1)*10./(TMath::Sqrt(sigN1_det*2*3.141593))*exp(-0.5*(TMath::Power((x[0]-x01),2.0)*sigN1(0,0)+TMath::Power(x[1]-y01,2.0)*sigN1(1,1)+2*sigN1(1,0)*(x[0]-x01)*(x[1]-y01)));

  double beamW1 = (1-nw_weight1)*10./(TMath::Sqrt(sigW1_det*2*3.141593))*exp(-0.5*(TMath::Power((x[0]-x01),2.0)*sigW1(0,0)+TMath::Power(x[1]-y01,2.0)*sigW1(1,1)+2*sigW1(1,0)*(x[0]-x01)*(x[1]-y01)));

  double beamN2 = (nw_weight2)*10./(TMath::Sqrt(sigN2_det*2*3.141593))*exp(-0.5*(TMath::Power((x[0]-x02),2.0)*sigN2(0,0)+TMath::Power(x[1]-y02,2.0)*sigN2(1,1)+2*sigN2(1,0)*(x[0]-x02)*(x[1]-y02)));

 double beamW2 = (1-nw_weight2)*10./(TMath::Sqrt(sigW2_det*2*3.141593))*exp(-0.5*(TMath::Power(x[0]-x02,2.0)*sigW2(0,0)+TMath::Power(x[1]-y02,2.0)*sigW2(1,1)+2*sigW2(1,0)*(x[0]-x02)*(x[1]-y02)));
  */

  double product = (beamN1+beamW1)*(beamN2+beamW2);

  double beamN1_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN1,2)))*TMath::Sqrt(TMath::Power(xwidthN1,2.))*TMath::Sqrt(TMath::Power(ywidthN1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN1,2))*(TMath::Power((x[0]-x01)/xwidthN1,2.0)+TMath::Power((x[1]-y01)/ywidthN1,2.0)-2*rhoN1*(x[0]-x01)*(x[1]-y01)/(xwidthN1*ywidthN1)));

  double beamW1_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW1,2)))*TMath::Sqrt(TMath::Power(xwidthW1,2.))*TMath::Sqrt(TMath::Power(ywidthW1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW1,2))*(TMath::Power((x[0]-x01)/xwidthW1,2.0)+TMath::Power((x[1]-y01)/ywidthW1,2.0)-2*rhoW1*(x[0]-x01)*(x[1]-y01)/(xwidthW1*ywidthW1)));

  double beamN2_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN2,2)))*TMath::Sqrt(TMath::Power(xwidthN2,2.))*TMath::Sqrt(TMath::Power(ywidthN2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN2,2))*(TMath::Power((x[0]-x02)/xwidthN2,2.0)+TMath::Power((x[1]-y02)/ywidthN2,2.0)-2*rhoN2*(x[0]-x02)*(x[1]-y02)/(xwidthN2*ywidthN2)));

  double beamW2_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW2,2)))*TMath::Sqrt(TMath::Power(xwidthW2,2.))*TMath::Sqrt(TMath::Power(ywidthW2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW2,2))*(TMath::Power((x[0]-x02)/xwidthW2,2.0)+TMath::Power((x[1]-y02)/ywidthW2,2.0)-2*rhoW2*(x[0]-x02)*(x[1]-y02)/(xwidthW2*ywidthW2)));

  double product2 = (nw_weight1*beamN1_ + (1.-nw_weight1)*beamW1_) * (nw_weight2*beamN2_ + (1.-nw_weight2)*beamW2_);

  return product2;

 
}


void vdmScanTreeAnalyzerDG_toys(TString suffix)
{
  TRandom3 r;
  r.SetSeed(0);
  int scansteps = 19.;
TFile *fAna = new TFile("VdmSimulatorFileTreeAnalyzedTOYS"+suffix+".root","recreate");
   TTree *outTrainTree = new TTree("outTrainTree","tree of VdM fit parameters");
  for(int z = 0 ; z < 10 ; z++){

  double in_yWidth1N = 1.6+r.Uniform(0.4);//1.6450//1.5987
  double in_xWidth1N = 1.6+r.Uniform(0.4);//1.7848//1.7249
  // double in_corr1N = -0.4+r.Uniform(0.8); //2.4951e-01//2.9898
  double in_corr1N = 0.;//-0.4+r.Uniform(0.8); //2.4951e-01//2.9898
  double in_yWidth1W = 2.0+r.Uniform(0.6);//2.4088//2.3419
  double in_xWidth1W = 2.0+r.Uniform(0.6);//2.3915//2.3524
  //double in_corr1W = -0.4+r.Uniform(0.8);//3.9284e-01//3.7407e-01
  double in_corr1W = 0.;//-0.4+r.Uniform(0.8);//3.9284e-01//3.7407e-01
  //double in_weight1 = r.Uniform(1.);//4.7091e-01//3.7466e-01
  double in_weight1 = 0.;//r.Uniform(1.);//4.7091e-01//3.7466e-01

  double in_yWidth2N = 1.6+r.Uniform(0.4);//1.6450//1.5987
  double in_xWidth2N = 1.6+r.Uniform(0.4);//1.7848//1.7249
  //double in_corr2N = -0.4+r.Uniform(0.8); //2.4951e-01//2.9898
  double in_corr2N = 0.;//-0.4+r.Uniform(0.8); //2.4951e-01//2.9898
  double in_yWidth2W = 2.0+r.Uniform(0.6);//2.4088//2.3419
  double in_xWidth2W = 2.0+r.Uniform(0.6);//2.3915//2.3524
  // double in_corr2W = -0.4+r.Uniform(0.8);//3.9284e-01//3.7407e-01
  double in_corr2W = 0.;//-0.4+r.Uniform(0.8);//3.9284e-01//3.7407e-01
  //double in_weight2 = r.Uniform(1.);//4.7091e-01//3.7466e-01
  double in_weight2 = 0.;// r.Uniform(1.);//4.7091e-01//3.7466e-01
 

  /* double in_yWidth2N =2.0;//1.9709//2.0000
  double in_xWidth2N = 1.8;//1.8197//1.8347
  double in_corr2N = -0.3;//-3.3854e-01//-3.9999
  double in_yWidth2W = 2.55;//2.4554//2.4361
  double in_xWidth2W = 2.3;//2.2129//2.1853
  double in_corr2W = -0.3;//-2.7845e-01//-2.5599e-01
  double in_weight2 = 0.4;//3.3235e-01//2.8638e-01*/
  cout<<in_corr2W<<endl;
  cout<<in_corr2N<<endl;
 cout<<in_corr1W<<endl;
 cout<<in_corr1N<<endl;
 cout<<in_weight2<<endl;
 cout<<in_weight1<<endl;
TF2 *multBeam = new TF2("multBeam",beamMultDG,-30,30,-30,30,18);





 multBeam->SetParameter(0, 0.0);
 multBeam->SetParameter(1, 0.0);
 multBeam->SetParameter(2,in_xWidth1N);
 multBeam->SetParameter(3,in_yWidth1N);
 multBeam->SetParameter(4,in_corr1N);
 multBeam->SetParameter(5,in_xWidth1W);
 multBeam->SetParameter(6,in_yWidth1W);
 multBeam->SetParameter(7,in_corr1W);
 multBeam->SetParameter(8,in_weight1);

 multBeam->SetParameter(9, 0.0);
 multBeam->SetParameter(10, 0.0);
 multBeam->SetParameter(11,in_xWidth2N);
 multBeam->SetParameter(12,in_yWidth2N);
 multBeam->SetParameter(13,in_corr2N);
 multBeam->SetParameter(14,in_xWidth2W);
 multBeam->SetParameter(15,in_yWidth2W);
 multBeam->SetParameter(16,in_corr2W);
 multBeam->SetParameter(17,in_weight2);
 
 cout<<"TEST INTEGRAL "<<multBeam->Integral(-30,30,-30,30)<<endl;

 double vertexResolution = 0.7;
 double nbinsxy = 40*scansteps;
 double beamSigma =2.;
 double beamSeperationMax = ((scansteps - 1.)/2.)*beamSigma/2.;
 double histLowEdge = (-beamSeperationMax-(beamSigma/2.));
 double histHighEdge = (beamSeperationMax+(beamSigma/2.));
 
 //TTree *treeVtxBeam1 = new TTree("treeVtxBeam1","tree of VdM vertices Beam1 at rest");
 //TTree *treeVtxBeam2 = new TTree("treeVtxBeam2","tree of VdM vertices Beam2 at rest");
 
 TH2F *Beam2MoveX_Add = new TH2F("Beam2MoveX_Add","Beam2MoveX_Add",nbinsxy,histLowEdge,histHighEdge, nbinsxy,histLowEdge, histHighEdge);
 TH2F *Beam2MoveY_Add = new TH2F("Beam2MoveY_Add","Beam2MoveY_Add",nbinsxy,histLowEdge,histHighEdge, nbinsxy,histLowEdge, histHighEdge);
 TH2F *Beam1MoveX_Add = new TH2F("Beam1MoveX_Add","Beam1MoveX_Add",nbinsxy,histLowEdge,histHighEdge, nbinsxy,histLowEdge, histHighEdge);
 TH2F *Beam1MoveY_Add = new TH2F("Beam1MoveY_Add","Beam1MoveY_Add",nbinsxy,histLowEdge,histHighEdge, nbinsxy,histLowEdge, histHighEdge);
 
 double vtxX = 0;
 double vtxY = 0;
 double vtxXsmeared = 0;
 double vtxYsmeared = 0;
 int scanXn = -100;
 int scanYn = -100;
 
 /*treeVtxBeam1->Branch("vtxX",&vtxX,"vtxX/D");
 treeVtxBeam1->Branch("vtxY",&vtxY,"vtxY/D");
 treeVtxBeam1->Branch("vtxXsmeared",&vtxXsmeared,"vtxXsmeared/D");
 treeVtxBeam1->Branch("vtxYsmeared",&vtxYsmeared,"vtxYsmeared/D");
 treeVtxBeam1->Branch("scanXn",&scanXn,"scanXn/I");
 treeVtxBeam1->Branch("scanYn",&scanYn,"scanYn/I");

 treeVtxBeam2->Branch("vtxX",&vtxX,"vtxX/D");
 treeVtxBeam2->Branch("vtxY",&vtxY,"vtxY/D");
 treeVtxBeam2->Branch("vtxXsmeared",&vtxXsmeared,"vtxXsmeared/D");
 treeVtxBeam2->Branch("vtxYsmeared",&vtxYsmeared,"vtxYsmeared/D");
 treeVtxBeam2->Branch("scanXn",&scanXn,"scanXn/I");
 treeVtxBeam2->Branch("scanYn",&scanYn,"scanYn/I");*/


 multBeam->SetNpy(500);
 multBeam->SetNpx(500);

 for(int scanX = 0; scanX<19; scanX+=1){
   scanXn = scanX;   
   multBeam->SetParameter(9,scanX-9);
   double integral = multBeam->Integral(-30,30,-30,30)*80;
   double integralP = r.PoissonD(integral);
   for(int n=0; n<integralP;n++){   
     multBeam->GetRandom2(vtxX,vtxY);
     vtxXsmeared = r.Gaus(vtxX,vertexResolution);
     vtxYsmeared = r.Gaus(vtxY,vertexResolution);
     Beam2MoveX_Add->Fill(vtxXsmeared,vtxYsmeared);

     //     treeVtxBeam1->Fill();
   }
 }

 multBeam->SetParameter(9,0.0);
 scanXn = -100;
 
 for(int scanY = 0; scanY<19; scanY+=1){
   scanYn = scanY;
   multBeam->SetParameter(10,scanY-9.);
   double integral = multBeam->Integral(-30,30,-30,30)*80;
   double integralP = r.PoissonD(integral);
   for(int n=0; n<integralP;n++){   
     multBeam->GetRandom2(vtxX,vtxY);
     vtxXsmeared = r.Gaus(vtxX,vertexResolution);
     vtxYsmeared = r.Gaus(vtxY,vertexResolution);
     Beam2MoveY_Add->Fill(vtxXsmeared,vtxYsmeared);
     //     treeVtxBeam1->Fill();
   }
 }

 scanYn = -100;
 multBeam->SetParameter(10,0.0);
 
 for(int scanX = 0; scanX<19; scanX+=1){
   scanXn = scanX;
   multBeam->SetParameter(0,scanX-9);
   double integral = multBeam->Integral(-30,30,-30,30)*80;
   double integralP = r.PoissonD(integral);
   for(int n=0; n<integralP;n++){
     multBeam->GetRandom2(vtxX,vtxY);
     vtxXsmeared = r.Gaus(vtxX,vertexResolution);
     vtxYsmeared = r.Gaus(vtxY,vertexResolution);
     Beam1MoveX_Add->Fill(vtxXsmeared,vtxYsmeared);
     //treeVtxBeam2->Fill();
   }
 }

 multBeam->SetParameter(0,0.0);
 scanXn = -100;
 
 for(int scanY = 0; scanY<19; scanY+=1){
   scanYn = scanY;
   multBeam->SetParameter(1,scanY-9.);
   double integral = multBeam->Integral(-30,30,-30,30)*80;
   double integralP = r.PoissonD(integral);
   for(int n=0; n<integralP;n++){    
     multBeam->GetRandom2(vtxX,vtxY);
     vtxXsmeared = r.Gaus(vtxX,vertexResolution);
     vtxYsmeared = r.Gaus(vtxY,vertexResolution); 
     Beam1MoveY_Add->Fill(vtxXsmeared,vtxYsmeared);
     //treeVtxBeam2->Fill();
   }
 }
 multBeam->SetParameter(1,0.0);


 std::cout<<Beam2MoveX_Add->GetEntries()<<std::endl;
	std::cout<<Beam2MoveY_Add->GetEntries()<<std::endl;
std::cout<<Beam1MoveX_Add->GetEntries()<<std::endl;
std::cout<<Beam1MoveY_Add->GetEntries()<<std::endl;
 



 //TCanvas *c1 = new TCanvas("c1","c1",400,600);
 // Beam2MoveX_Add->Draw();


 RooRealVar xVar("xVar","xVar",-10.0,10.0) ;
 RooRealVar yVar("yVar","yVar",-10.0,10.0) ;
 xVar.setBins(10000,"cache");
 yVar.setBins(10000,"cache");
            
 RooRealVar yWidthN1("yWidthN1","yWidthN1",0.7,1.89) ;
 RooRealVar xWidthN1("xWidthN1","xWidthN1",0.7,1.89) ;  
 RooRealVar rho_N1("rho_N1","rho_N1",-0.49,0.49) ;
 RooRealVar yWidthW1("yWidthW1","yWidthW1",1.9,3.6) ;
 RooRealVar xWidthW1("xWidthW1","xWidthW1",1.9,3.6) ;  
 RooRealVar rho_W1("rho_W1","rho_W1",-0.49,0.49) ;
 RooRealVar w1("w1","w1",0.0,1.0) ;

 RooRealVar yWidthN2("yWidthN2","yWidthN2",0.7,1.89) ;
 RooRealVar xWidthN2("xWidthN2","xWidthN2",0.7,1.89) ;  
 RooRealVar rho_N2("rho_N2","rho_N2",-0.49,0.49) ;
 RooRealVar yWidthW2("yWidthW2","yWidthW2",1.9,3.6) ;
 RooRealVar xWidthW2("xWidthW2","xWidthW2",1.9,3.6) ;  
 RooRealVar rho_W2("rho_W2","rho_W2",-0.49,0.49) ;
 RooRealVar w2("w2","w2",0.0,1.0) ;

 MyPdfV3  beam1RestVerticesUnfold_XScan("beam1RestVerticesUnfold_Xscan","beam1RestVerticesUnfold_Xscan",xVar,yVar,w1,rho_N1,xWidthN1,yWidthN1,rho_W1,xWidthW1,yWidthW1,w2,yWidthN2,yWidthW2);
 MyPdfV4  beam1RestVerticesUnfold_YScan("beam1RestVerticesUnfold_Yscan","beam1RestVerticesUnfold_Yscan",xVar,yVar,w1,rho_N1,xWidthN1,yWidthN1,rho_W1,xWidthW1,yWidthW1,w2,xWidthN2,xWidthW2);
 MyPdfV3  beam2RestVerticesUnfold_XScan("beam2RestVerticesUnfold_Xscan","beam2RestVerticesUnfold_Xscan",xVar,yVar,w2,rho_N2,xWidthN2,yWidthN2,rho_W2,xWidthW2,yWidthW2,w1,yWidthN1,yWidthW1);
 MyPdfV4  beam2RestVerticesUnfold_YScan("beam2RestVerticesUnfold_Yscan","beam2RestVerticesUnfold_Yscan",xVar,yVar,w2,rho_N2,xWidthN2,yWidthN2,rho_W2,xWidthW2,yWidthW2,w1,xWidthN1,xWidthW1);

 RooDataHist *scanXBeam1RestDataHist = new RooDataHist("scanXBeam1RestDataHist","scanXBeam1RestDataHist",RooArgList(xVar,yVar),Beam2MoveX_Add);
 RooDataHist *scanYBeam1RestDataHist = new RooDataHist("scanYBeam1RestDataHist","scanYBeam1RestDataHist",RooArgList(xVar,yVar),Beam2MoveY_Add);
 RooDataHist *scanXBeam2RestDataHist = new RooDataHist("scanXBeam2RestDataHist","scanXBeam2RestDataHist",RooArgList(xVar,yVar),Beam1MoveX_Add);
 RooDataHist *scanYBeam2RestDataHist = new RooDataHist("scanYBeam2RestDataHist","scanYBeam2RestDataHist",RooArgList(xVar,yVar),Beam1MoveY_Add);

 RooCategory sample("sample","sample") ;
 sample.defineType("X_ScanData_Beam1Rest") ;
 sample.defineType("Y_ScanData_Beam1Rest") ;
 sample.defineType("X_ScanData_Beam2Rest") ;
 sample.defineType("Y_ScanData_Beam2Rest") ;
 RooDataHist combData("combData","combined data",RooArgSet(xVar,yVar),Index(sample),Import("X_ScanData_Beam1Rest",*scanXBeam1RestDataHist),Import("Y_ScanData_Beam1Rest",*scanYBeam1RestDataHist),Import("X_ScanData_Beam2Rest",*scanXBeam2RestDataHist),Import("Y_ScanData_Beam2Rest",*scanYBeam2RestDataHist)); 
 RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
 simPdf.addPdf(beam1RestVerticesUnfold_XScan,"X_ScanData_Beam1Rest");
 simPdf.addPdf(beam1RestVerticesUnfold_YScan,"Y_ScanData_Beam1Rest");
 simPdf.addPdf(beam2RestVerticesUnfold_XScan,"X_ScanData_Beam2Rest");
 simPdf.addPdf(beam2RestVerticesUnfold_YScan,"Y_ScanData_Beam2Rest");
 RooFitResult *result = simPdf.fitTo(combData,RooFit::NumCPU(5),RooFit::Save()/*,RooFit::PrintLevel(3),RooFit::Verbose(1)*/);


 std::cout<<"BLUB"<<in_yWidth2N<<endl;
 // result->Print();

 /*TFile *fAna = new TFile("VdmSimulatorFileTreeAnalyzedTOYS"+suffix+".root","recreate");
   TTree *outTrainTree = new TTree("outTrainTree","tree of VdM fit parameters");*/

 double yWidth1N_true = in_yWidth1N;
 double xWidth1N_true = in_xWidth1N;
 double corr1N_true = in_corr1N;
 double yWidth1W_true = in_yWidth1W;
 double xWidth1W_true = in_xWidth1W;
 double corr1W_true = in_corr1W;
 double weight1_true = in_weight1;
 
 double yWidth2N_true = in_yWidth2N;
 double xWidth2N_true = in_xWidth2N;
 double corr2N_true = in_corr2N;
 double yWidth2W_true = in_yWidth2W;
 double xWidth2W_true = in_xWidth2W;
 double corr2W_true = in_corr2W;
 double weight2_true = in_weight2;
                       
 double yWidth1N_fit = yWidthN1.getValV();
 double xWidth1N_fit = xWidthN1.getValV();
 double corr1N_fit = rho_N1.getValV();
 double yWidth1W_fit = yWidthW1.getValV();
 double xWidth1W_fit = xWidthW1.getValV();
 double corr1W_fit = rho_W1.getValV();
 double weight1_fit = w1.getValV();
 
 double yWidth2N_fit = yWidthN2.getValV();
 double xWidth2N_fit = xWidthN2.getValV();
 double corr2N_fit = rho_N2.getValV();
 double yWidth2W_fit = yWidthW2.getValV();
 double xWidth2W_fit = xWidthW2.getValV();
 double corr2W_fit = rho_W2.getValV();
 double weight2_fit = w2.getValV();
 
 double yWidth1N_fitErr = yWidthN1.getError();
 double xWidth1N_fitErr = xWidthN1.getError();
 double corr1N_fitErr = rho_N1.getError();
 double yWidth1W_fitErr = yWidthW1.getError();
 double xWidth1W_fitErr = xWidthW1.getError();
 double corr1W_fitErr = rho_W1.getError();
 double weight1_fitErr = w1.getError();
 
 double yWidth2N_fitErr = yWidthN2.getError();
 double xWidth2N_fitErr = xWidthN2.getError();
 double corr2N_fitErr = rho_N2.getError();
 double yWidth2W_fitErr = yWidthW2.getError();
 double xWidth2W_fitErr = xWidthW2.getError();
 double corr2W_fitErr = rho_W2.getError();
 double weight2_fitErr = w2.getError();

 double yWidth1N_diff = (yWidth1N_fit - yWidth1N_true)/yWidth1N_true;
 double xWidth1N_diff = (xWidth1N_fit - xWidth1N_true)/xWidth1N_true;
 double corr1N_diff = (corr1N_fit - corr1N_true)/TMath::Abs(corr1N_true);
 double yWidth1W_diff = (yWidth1W_fit - yWidth1W_true)/yWidth1W_true;
 double xWidth1W_diff = (xWidth1W_fit - xWidth1W_true)/xWidth1W_true;
 double corr1W_diff = (corr1W_fit - corr1W_true)/TMath::Abs(corr1W_true);
 double weight1_diff = (weight1_fit - weight1_true)/weight1_true;

 double yWidth2N_diff = (yWidth2N_fit - yWidth2N_true)/yWidth2N_true;
 double xWidth2N_diff = (xWidth2N_fit - xWidth2N_true)/xWidth2N_true;
 double corr2N_diff = (corr2N_fit - corr2N_true)/TMath::Abs(corr2N_true);
 double yWidth2W_diff = (yWidth2W_fit - yWidth2W_true)/yWidth2W_true;
 double xWidth2W_diff = (xWidth2W_fit - xWidth2W_true)/xWidth2W_true;
 double corr2W_diff = (corr2W_fit - corr2W_true)/TMath::Abs(corr2W_true);
 double weight2_diff = (weight2_fit - weight2_true)/weight2_true;
 
 TF2 *multBeamT = new TF2("multBeamT",beamMultDGt,-30,30,-30,30,18);
 
 multBeamT->SetParameter(0, 0.0);
 multBeamT->SetParameter(1, 0.0);
 multBeamT->SetParameter(2,in_xWidth1N);
 multBeamT->SetParameter(3,in_yWidth1N);
 multBeamT->SetParameter(4,in_corr1N);
 multBeamT->SetParameter(5,in_xWidth1W);
 multBeamT->SetParameter(6,in_yWidth1W);
 multBeamT->SetParameter(7,in_corr1W);
 multBeamT->SetParameter(8,in_weight1);
 multBeamT->SetParameter(9, 0.0);
 multBeamT->SetParameter(10, 0.0);
 multBeamT->SetParameter(11,in_xWidth2N);
 multBeamT->SetParameter(12,in_yWidth2N);
 multBeamT->SetParameter(13,in_corr2N);
 multBeamT->SetParameter(14,in_xWidth2W);
 multBeamT->SetParameter(15,in_yWidth2W);
 multBeamT->SetParameter(16,in_corr2W);
 multBeamT->SetParameter(17,in_weight2);

 double overlapTrue = multBeamT->Integral(-30,30,-30,30);

 multBeamT->SetParameter(0, 0.0);
 multBeamT->SetParameter(1, 0.0);
 multBeamT->SetParameter(2,xWidth1N_fit);
 multBeamT->SetParameter(3,yWidth1N_fit);
 multBeamT->SetParameter(4,corr1N_fit);
 multBeamT->SetParameter(5,xWidth1W_fit);
 multBeamT->SetParameter(6,yWidth1W_fit);
 multBeamT->SetParameter(7,corr1W_fit);
 multBeamT->SetParameter(8,weight1_fit);
 multBeamT->SetParameter(9, 0.0);
 multBeamT->SetParameter(10, 0.0);
 multBeamT->SetParameter(11,xWidth2N_fit);
 multBeamT->SetParameter(12,yWidth2N_fit);
 multBeamT->SetParameter(13,corr2N_fit);
 multBeamT->SetParameter(14,xWidth2W_fit);
 multBeamT->SetParameter(15,yWidth2W_fit);
 multBeamT->SetParameter(16,corr2W_fit);
 multBeamT->SetParameter(17,weight2_fit);

 double overlapFit = multBeamT->Integral(-30,30,-30,30);
 double overlapDiff = (overlapFit - overlapTrue)/overlapTrue;

 outTrainTree->Branch("overlapTrue",&overlapTrue,"overlapTrue/D");
 outTrainTree->Branch("overlapFit",&overlapFit,"overlapFit/D");
 outTrainTree->Branch("overlapDiff",&overlapDiff,"overlapDiff/D");

 outTrainTree->Branch("yWidth1N_true",&yWidth1N_true,"yWidth1N_true/D");
 outTrainTree->Branch("xWidth1N_true",&xWidth1N_true,"xWidth1N_true/D");
 outTrainTree->Branch("corr1N_true",&corr1N_true,"corr1N_true/D");
 outTrainTree->Branch("yWidth1W_true",&yWidth1W_true,"yWidth1W_true/D");
 outTrainTree->Branch("xWidth1W_true",&xWidth1W_true,"xWidth1W_true/D");
 outTrainTree->Branch("corr1W_true",&corr1W_true,"corr1W_true/D");
 outTrainTree->Branch("weight1_true",&weight1_true,"weight1_true/D");

 outTrainTree->Branch("yWidth2N_true",&yWidth2N_true,"yWidth2N_true/D");
 outTrainTree->Branch("xWidth2N_true",&xWidth2N_true,"xWidth2N_true/D");
 outTrainTree->Branch("corr2N_true",&corr2N_true,"corr2N_true/D");
 outTrainTree->Branch("yWidth2W_true",&yWidth2W_true,"yWidth2W_true/D");
 outTrainTree->Branch("xWidth2W_true",&xWidth2W_true,"xWidth2W_true/D");
 outTrainTree->Branch("corr2W_true",&corr2W_true,"corr2W_true/D");
 outTrainTree->Branch("weight2_true",&weight2_true,"weight2_true/D");

 outTrainTree->Branch("yWidth1N_fit",&yWidth1N_fit,"yWidth1N_fit/D");
 outTrainTree->Branch("xWidth1N_fit",&xWidth1N_fit,"xWidth1N_fit/D");
 outTrainTree->Branch("corr1N_fit",&corr1N_fit,"corr1N_fit/D");
 outTrainTree->Branch("yWidth1W_fit",&yWidth1W_fit,"yWidth1W_fit/D");
 outTrainTree->Branch("xWidth1W_fit",&xWidth1W_fit,"xWidth1W_fit/D");
 outTrainTree->Branch("corr1W_fit",&corr1W_fit,"corr1W_fit/D");
 outTrainTree->Branch("weight1_fit",&weight1_fit,"weight1_fit/D");

 outTrainTree->Branch("yWidth2N_fit",&yWidth2N_fit,"yWidth2N_fit/D");
 outTrainTree->Branch("xWidth2N_fit",&xWidth2N_fit,"xWidth2N_fit/D");
 outTrainTree->Branch("corr2N_fit",&corr2N_fit,"corr2N_fit/D");
 outTrainTree->Branch("yWidth2W_fit",&yWidth2W_fit,"yWidth2W_fit/D");
 outTrainTree->Branch("xWidth2W_fit",&xWidth2W_fit,"xWidth2W_fit/D");
 outTrainTree->Branch("corr2W_fit",&corr2W_fit,"corr2W_fit/D");
 outTrainTree->Branch("weight2_fit",&weight2_fit,"weight2_fit/D");



 outTrainTree->Branch("yWidth1N_fitErr",&yWidth1N_fitErr,"yWidth1N_fitErr/D");
 outTrainTree->Branch("xWidth1N_fitErr",&xWidth1N_fitErr,"xWidth1N_fitErr/D");
 outTrainTree->Branch("corr1N_fitErr",&corr1N_fitErr,"corr1N_fitErr/D");
 outTrainTree->Branch("yWidth1W_fitErr",&yWidth1W_fitErr,"yWidth1W_fitErr/D");
 outTrainTree->Branch("xWidth1W_fitErr",&xWidth1W_fitErr,"xWidth1W_fitErr/D");
 outTrainTree->Branch("corr1W_fitErr",&corr1W_fitErr,"corr1W_fitErr/D");
 outTrainTree->Branch("weight1_fitErr",&weight1_fitErr,"weight1_fitErr/D");

 outTrainTree->Branch("yWidth2N_fitErr",&yWidth2N_fitErr,"yWidth2N_fitErr/D");
 outTrainTree->Branch("xWidth2N_fitErr",&xWidth2N_fitErr,"xWidth2N_fitErr/D");
 outTrainTree->Branch("corr2N_fitErr",&corr2N_fitErr,"corr2N_fitErr/D");
 outTrainTree->Branch("yWidth2W_fitErr",&yWidth2W_fitErr,"yWidth2W_fitErr/D");
 outTrainTree->Branch("xWidth2W_fitErr",&xWidth2W_fitErr,"xWidth2W_fitErr/D");
 outTrainTree->Branch("corr2W_fitErr",&corr2W_fitErr,"corr2W_fitErr/D");
 outTrainTree->Branch("weight2_fitErr",&weight2_fitErr,"weight2_fitErr/D");


 outTrainTree->Branch("yWidth1N_diff",&yWidth1N_diff,"yWidth1N_diff/D");
 outTrainTree->Branch("xWidth1N_diff",&xWidth1N_diff,"xWidth1N_diff/D");
 outTrainTree->Branch("corr1N_diff",&corr1N_diff,"corr1N_diff/D");
 outTrainTree->Branch("yWidth1W_diff",&yWidth1W_diff,"yWidth1W_diff/D");
 outTrainTree->Branch("xWidth1W_diff",&xWidth1W_diff,"xWidth1W_diff/D");
 outTrainTree->Branch("corr1W_diff",&corr1W_diff,"corr1W_diff/D");
 outTrainTree->Branch("weight1_diff",&weight1_diff,"weight1_diff/D");

 outTrainTree->Branch("yWidth2N_diff",&yWidth2N_diff,"yWidth2N_diff/D");
 outTrainTree->Branch("xWidth2N_diff",&xWidth2N_diff,"xWidth2N_diff/D");
 outTrainTree->Branch("corr2N_diff",&corr2N_diff,"corr2N_diff/D");
 outTrainTree->Branch("yWidth2W_diff",&yWidth2W_diff,"yWidth2W_diff/D");
 outTrainTree->Branch("xWidth2W_diff",&xWidth2W_diff,"xWidth2W_diff/D");
 outTrainTree->Branch("corr2W_diff",&corr2W_diff,"corr2W_diff/D");
 outTrainTree->Branch("weight2_diff",&weight2_diff,"weight2_diff/D");
 

 //for


 outTrainTree->Fill();
  }
 outTrainTree->Write();
 
															   

}
