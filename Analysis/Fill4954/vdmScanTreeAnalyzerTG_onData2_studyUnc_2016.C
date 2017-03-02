//#include "MyPdfV1.h"
#include "../../beamdensities/TripleGauss_V1.h"
#include "../../beamdensities/TripleGauss_V2.h"
//#include "gauss2dProjectionV4.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TMath.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <sstream>
#include "TRint.h"
#include "TTree.h"
#include "TObject.h"
//#include "TVirtualFFT.h"
//#include "RooGlobalFunc.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "TFormula.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction2Binding.h"
#include "RooTFnBinding.h"
#include "RooMultiVarGaussian.h"
#include "RooProdPdf.h"
//#include "Roo2DimGauss.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
//#include "RooWorkspace.h"
#include "TStyle.h"
//#include "TMVA/Tools.h"
#include "../../TMVA-v4.2.0/TMVA/Reader.h"
#include "../../TMVA-v4.2.0/TMVA/Tools.h"
using namespace RooFit;


Double_t beamMultTG(Double_t *x,Double_t *par)
{
  Double_t arg = 0;
  Double_t pi = 3.1415926;

  double x01 = par[0];
  double y01 = par[1];
  double xwidthN1 = par[2];
  double ywidthN1 = par[3];
  double rhoN1 = par[4];
  double xwidthM1 = par[5];
  double ywidthM1 = par[6];
  double rhoM1 = par[7];
  double xwidthW1 = par[8];
  double ywidthW1 = par[9];
  double rhoW1 = par[10];
  double nw_weight1N = par[11];
  double nw_weight1M = par[12];
  double nw_weight1W = 1 - par[11] - par[12];

  double x02 = par[13];
  double y02 = par[14];
  double xwidthN2 = par[15];
  double ywidthN2 = par[16];
  double rhoN2 = par[17];
  double xwidthM2 = par[18];
  double ywidthM2 = par[19];
  double rhoM2 = par[20];
  double xwidthW2 = par[21];
  double ywidthW2 = par[22];
  double rhoW2 = par[23];
  double nw_weight2N = par[24];
  double nw_weight2M = par[25];
  double nw_weight2W = 1 - par[24] - par[25];

  double beamN1 = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN1,2)))*TMath::Sqrt(TMath::Power(xwidthN1,2.))*TMath::Sqrt(TMath::Power(ywidthN1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN1,2))*(TMath::Power((x[0]-x01)/xwidthN1,2.0)+TMath::Power((x[1]-y01)/ywidthN1,2.0)-2*rhoN1*(x[0]-x01)*(x[1]-y01)/(xwidthN1*ywidthN1)));
  double beamM1 = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoM1,2)))*TMath::Sqrt(TMath::Power(xwidthM1,2.))*TMath::Sqrt(TMath::Power(ywidthM1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoM1,2))*(TMath::Power((x[0]-x01)/xwidthM1,2.0)+TMath::Power((x[1]-y01)/ywidthM1,2.0)-2*rhoM1*(x[0]-x01)*(x[1]-y01)/(xwidthM1*ywidthM1)));
  double beamW1 = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW1,2)))*TMath::Sqrt(TMath::Power(xwidthW1,2.))*TMath::Sqrt(TMath::Power(ywidthW1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW1,2))*(TMath::Power((x[0]-x01)/xwidthW1,2.0)+TMath::Power((x[1]-y01)/ywidthW1,2.0)-2*rhoW1*(x[0]-x01)*(x[1]-y01)/(xwidthW1*ywidthW1)));
  double beam1 = nw_weight1N * beamN1 + nw_weight1M * beamM1 + nw_weight1W * beamW1;

  double beamN2 = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN2,2)))*TMath::Sqrt(TMath::Power(xwidthN2,2.))*TMath::Sqrt(TMath::Power(ywidthN2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN2,2))*(TMath::Power((x[0]-x02)/xwidthN2,2.0)+TMath::Power((x[1]-y02)/ywidthN2,2.0)-2*rhoN2*(x[0]-x02)*(x[1]-y02)/(xwidthN2*ywidthN2)));
  double beamM2 = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoM2,2)))*TMath::Sqrt(TMath::Power(xwidthM2,2.))*TMath::Sqrt(TMath::Power(ywidthM2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoM2,2))*(TMath::Power((x[0]-x02)/xwidthM2,2.0)+TMath::Power((x[1]-y02)/ywidthM2,2.0)-2*rhoM2*(x[0]-x02)*(x[1]-y02)/(xwidthM2*ywidthM2)));
  double beamW2 = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW2,2)))*TMath::Sqrt(TMath::Power(xwidthW2,2.))*TMath::Sqrt(TMath::Power(ywidthW2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW2,2))*(TMath::Power((x[0]-x02)/xwidthW2,2.0)+TMath::Power((x[1]-y02)/ywidthW2,2.0)-2*rhoW2*(x[0]-x02)*(x[1]-y02)/(xwidthW2*ywidthW2)));
  double beam2 = nw_weight2N * beamN2 + nw_weight2M * beamM2 + nw_weight2W * beamW2;

  double product = beam1 * beam2;

  return product;


}


void vdmScanTreeAnalyzerTG_onData2_studyUnc_2016(TString bx)
{

    double scaling = 0.00458;

  //TFile *f = TFile::Open("testNEW8_ext.root");
  //TFile *f = TFile::Open("../vdm_ReRecoAnalysis/newBeamImaging.root");
//TFile *f = TFile::Open("2016Scans_v5.root");



  //TFile *f = TFile::Open("../vdm_ReRecoAnalysis/newbemimagingfine.root");

   TFile *f = TFile::Open("/eos/cms/store/user/jsalfeld/vdmScan_2016ReRecoJan/mergedJan172016.root");

   TString suff="StronRescale";

//TFile *f = TFile::Open("testNEW8_coarse_ext.root");
   //TString bunchStr[5] = {"51","771","1631","2211","2674"};
//TString bunchStr[5] = {"41","281","872","1783","2063"};
TString bunchStr[1] = {bx};
 gStyle->SetOptStat(0);


/*TMVA::Tools::Instance();

TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );*>



 float yWidth1N_fitf;//1.6450//1.5987
  float xWidth1N_fitf ;//1.7848//1.7249
  float corr1N_fitf ; //2.4951e-01//2.9898
  float yWidth1W_fitf ;//2.4088//2.3419
  float xWidth1W_fitf ;//2.3915//2.3524
  float corr1W_fitf ;//3.9284e-01//3.7407e-01
  float weight1_fitf;//4.7091e-01//3.7466e-01

  float yWidth2N_fitf ;//1.6450//1.5987
  float xWidth2N_fitf ;//1.7848//1.7249
  float corr2N_fitf ; //2.4951e-01//2.9898
  float yWidth2W_fitf ;//2.4088//2.3419
  float xWidth2W_fitf ;//2.3915//2.3524
  float corr2W_fitf ;//3.9284e-01//3.7407e-01
  float weight2_fitf ;//


  reader->AddVariable("yWidth1N_fit",&yWidth1N_fitf );
  reader->AddVariable("xWidth1N_fit",&xWidth1N_fitf );
  reader->AddVariable("corr1N_fit",&corr1N_fitf );
  reader->AddVariable("yWidth1W_fit",&yWidth1W_fitf );
  reader->AddVariable("xWidth1W_fit",&xWidth1W_fitf );
  reader->AddVariable("corr1W_fit",&corr1W_fitf );
  reader->AddVariable("weight1_fit",&weight1_fitf );

  reader->AddVariable("yWidth2N_fit",&yWidth2N_fitf );
  reader->AddVariable("xWidth2N_fit",&xWidth2N_fitf );
  reader->AddVariable("corr2N_fit",&corr2N_fitf );
  reader->AddVariable("yWidth2W_fit",&yWidth2W_fitf );
  reader->AddVariable("xWidth2W_fit",&xWidth2W_fitf );
  reader->AddVariable("corr2W_fit",&corr2W_fitf );
  reader->AddVariable("weight2_fit",&weight2_fitf );

  reader->BookMVA( "MLP", "../../TMVA-v4.2.0/test/weights/TMVARegression_MLP.weights.xml" );
*/

int i = 0;
 //for(int i=0;i<5;i++){

 TH2F *Beam2MoveX_Add = (TH2F*) f->Get("Beam2MoveX_bunch"+bunchStr[i]+"Add");
 TH2F *Beam2MoveY_Add = (TH2F*) f->Get("Beam2MoveY_bunch"+bunchStr[i]+"Add");
 TH2F *Beam1MoveX_Add = (TH2F*) f->Get("Beam1MoveX_bunch"+bunchStr[i]+"Add");
 TH2F *Beam1MoveY_Add = (TH2F*) f->Get("Beam1MoveY_bunch"+bunchStr[i]+"Add");



 cout<<endl<<endl<<endl<<"*******************"<<"Processing"<<Beam2MoveX_Add->GetTitle()<<"*****************"<<endl<<endl<<endl;



 RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-7) ;
 RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-7) ;

 RooRealVar y0_1("y0_1","y0_1",-0.5,0.5) ;
 RooRealVar x0_1("x0_1","x0_1",-0.5,0.5) ;
 RooRealVar y0_2("y0_2","y0_2",-0.5,0.5) ;
 RooRealVar x0_2("x0_2","x0_2",-0.5,0.5) ;

 RooRealVar y0_12("y0_12","y0_12",-0.5,0.5) ;
 RooRealVar x0_12("x0_12","x0_12",-0.5,0.5) ;
 RooRealVar y0_22("y0_22","y0_22",-0.5,0.5) ;
 RooRealVar x0_22("x0_22","x0_22",-0.5,0.5) ;

  RooRealVar xVar("xVar","xVar",-10.0,10.0) ;
  RooRealVar yVar("yVar","yVar",-10.0,10.0) ;
  xVar.setBins(10000,"cache");
  yVar.setBins(10000,"cache");



  /* RooRealVar yWidthN1("yWidthN1","yWidthN1",1.5,2.0) ;
  RooRealVar xWidthN1("xWidthN1","xWidthN1",1.5,2.0) ;
  RooRealVar rho_N1("rho_N1","rho_N1",-0.42,0.42) ;
  RooRealVar yWidthW1("yWidthW1","yWidthW1",1.8,2.7) ;
  RooRealVar xWidthW1("xWidthW1","xWidthW1",1.8,2.7) ;
  RooRealVar rho_W1("rho_W1","rho_W1",-0.42,0.42) ;
  RooRealVar w1("w1","w1",-0.5,1.0) ;

  RooRealVar yWidthN2("yWidthN2","yWidthN2",1.5,2.0) ;
  RooRealVar xWidthN2("xWidthN2","xWidthN2",1.5,2.0) ;
  RooRealVar rho_N2("rho_N2","rho_N2",-0.4,0.45) ;
  RooRealVar yWidthW2("yWidthW2","yWidthW2",1.8,2.7) ;
  RooRealVar xWidthW2("xWidthW2","xWidthW2",1.8,2.7) ;
  RooRealVar rho_W2("rho_W2","rho_W2",-0.32,0.3) ;
  RooRealVar w2("w2","w2",-0.5,1.0) ;*/






  /*RooRealVar yWidthN1("yWidthN1","yWidthN1",1.1,2.7) ;
  RooRealVar xWidthN1("xWidthN1","xWidthN1",1.1,2.7) ;
  RooRealVar rho_N1("rho_N1","rho_N1",-0.42,0.42) ;
  RooRealVar yWidthW1("yWidthW1","yWidthW1",1.5,8.7) ;
  RooRealVar xWidthW1("xWidthW1","xWidthW1",1.5,8.7) ;
  RooRealVar rho_W1("rho_W1","rho_W1",-0.45,0.45) ;
  RooRealVar w1("w1","w1",0.0,1.0) ;

  RooRealVar yWidthN2("yWidthN2","yWidthN2",1.1,2.7) ;
  RooRealVar xWidthN2("xWidthN2","xWidthN2",1.1,2.7) ;
  RooRealVar rho_N2("rho_N2","rho_N2",-0.4,0.45) ;
  RooRealVar yWidthW2("yWidthW2","yWidthW2",1.5,8.7) ;
  RooRealVar xWidthW2("xWidthW2","xWidthW2",1.5,8.7) ;
  RooRealVar rho_W2("rho_W2","rho_W2",-0.45,0.45) ;
  RooRealVar w2("w2","w2",0.0,1.0) ;
  */

  double pi = 3.141592654;

  RooRealVar yWidthN1("yWidthN1","yWidthN1",1.3,3.0) ;
  RooRealVar xWidthN1("xWidthN1","xWidthN1",1.3,3.0) ;
  RooRealVar rho_N1("rho_N1","rho_N1",-0.48,0.48) ;
  //RooRealVar yWidthM1("yWidthM1","yWidthM1",1.3,3.0) ;
  //RooRealVar xWidthM1("xWidthM1","xWidthM1",1.3,3.0) ;
  RooRealVar yWidthM1Diff("yWidthM1Diff","yWidthM1Diff",0.01,1.2);
  RooRealVar xWidthM1Diff("xWidthM1Diff","xWidthM1Diff",0.01,1.2);
  RooFormulaVar yWidthM1("yWidthM1","yWidthN1+yWidthM1Diff",RooArgSet(yWidthN1,yWidthM1Diff));
  RooFormulaVar xWidthM1("xWidthM1","xWidthN1+xWidthM1Diff",RooArgSet(xWidthN1,xWidthM1Diff));
  RooRealVar rho_M1("rho_M1","rho_M1",-0.48,0.48) ;
  //RooRealVar yWidthW1("yWidthW1","yWidthW1",1.3,3.0) ;
  //RooRealVar xWidthW1("xWidthW1","xWidthW1",1.3,3.0) ;
  RooRealVar yWidthW1Diff("yWidthW1Diff","yWidthW1Diff",0.01,1.2);
  RooRealVar xWidthW1Diff("xWidthW1Diff","xWidthW1Diff",0.01,1.2);
  RooFormulaVar yWidthW1("yWidthW1","yWidthN1+yWidthM1Diff+yWidthW1Diff",RooArgSet(yWidthN1,yWidthM1Diff,yWidthW1Diff));
  RooFormulaVar xWidthW1("xWidthW1","xWidthN1+xWidthM1Diff+xWidthW1Diff",RooArgSet(xWidthN1,xWidthM1Diff,xWidthW1Diff));
  RooRealVar rho_W1("rho_W1","rho_W1",-0.48,0.48) ;
  RooRealVar theta1("theta1","theta1",0.0,0.5*pi) ;
  RooRealVar phi1("phi1","phi1",0.0,0.5*pi) ;
  RooFormulaVar w1N("w1N","sin(theta1)**2*cos(phi1)**2", RooArgSet(theta1,phi1)) ;
  RooFormulaVar w1M("w1M","sin(theta1)**2*sin(phi1)**2", RooArgSet(theta1,phi1)) ;
  RooFormulaVar w1W("w1W","cos(theta1)**2", RooArgSet(theta1)) ;

  RooRealVar yWidthN2("yWidthN2","yWidthN2",1.3,3.0) ;
  RooRealVar xWidthN2("xWidthN2","xWidthN2",1.3,3.0) ;
  RooRealVar rho_N2("rho_N2","rho_N2",-0.48,0.48) ;
  //RooRealVar yWidthM2("yWidthM2","yWidthM2",1.3,3.0) ;
  //RooRealVar xWidthM2("xWidthM2","xWidthM2",1.3,3.0) ;
  RooRealVar yWidthM2Diff("yWidthM2Diff","yWidthM2Diff",0.01,1.2);
  RooRealVar xWidthM2Diff("xWidthM2Diff","xWidthM2Diff",0.01,1.2);
  RooFormulaVar yWidthM2("yWidthM2","yWidthN2+yWidthM2Diff",RooArgSet(yWidthN2,yWidthM2Diff));
  RooFormulaVar xWidthM2("xWidthM2","xWidthN2+xWidthM2Diff",RooArgSet(xWidthN2,xWidthM2Diff));
  RooRealVar rho_M2("rho_M2","rho_M2",-0.48,0.48) ;
  //RooRealVar yWidthW2("yWidthW2","yWidthW2",1.3,3.0) ;
  //RooRealVar xWidthW2("xWidthW2","xWidthW2",1.3,3.0) ;
  RooRealVar yWidthW2Diff("yWidthW2Diff","yWidthW2Diff",0.01,1.2);
  RooRealVar xWidthW2Diff("xWidthW2Diff","xWidthW2Diff",0.01,1.2);
  RooFormulaVar yWidthW2("yWidthW2","yWidthN2+yWidthM2Diff+yWidthW2Diff",RooArgSet(yWidthN2,yWidthM2Diff,yWidthW2Diff));
  RooFormulaVar xWidthW2("xWidthW2","xWidthN2+xWidthM2Diff+xWidthW2Diff",RooArgSet(xWidthN2,xWidthM2Diff,xWidthW2Diff));
  RooRealVar rho_W2("rho_W2","rho_W2",-0.48,0.48) ;
  RooRealVar theta2("theta2","theta2",0.0,0.5*pi) ;
  RooRealVar phi2("phi2","phi2",0.0,0.5*pi) ;
  RooFormulaVar w2N("w2N","sin(theta2)**2*cos(phi2)**2", RooArgSet(theta2,phi2)) ;
  RooFormulaVar w2M("w2M","sin(theta2)**2*sin(phi2)**2", RooArgSet(theta2,phi2)) ;
  RooFormulaVar w2W("w2W","cos(theta2)**2", RooArgSet(theta2)) ;

  RooRealVar vtxRes("vtxRes","vtxRes",0.00356/scaling) ;
  vtxRes.setConstant();

  RooRealVar xWidthRes("xWidthRes","xWidthRes",0.62,0.62) ;
  RooRealVar yWidthRes("yWidthRes","yWidthRes",0.62,0.62) ;
  RooRealVar meanRes("meanRes","meanRes",0.0,0.0) ;



  RooGaussian resX("resX","resX",xVar,meanRes,xWidthRes);
  RooGaussian resY("resY","resY",yVar,meanRes,yWidthRes);


  TripleGauss_V1  beam1RestVerticesUnfold_XScan("beam1RestVerticesUnfold_Xscan","beam1RestVerticesUnfold_Xscan",xVar,yVar,x0_1,y0_1,w1N,w1M,rho_N1,xWidthN1,yWidthN1,rho_M1,xWidthM1,yWidthM1,rho_W1,xWidthW1,yWidthW1,w2N,w2M,yWidthN2,xWidthM2,yWidthW2,vtxRes);
  TripleGauss_V2  beam1RestVerticesUnfold_YScan("beam1RestVerticesUnfold_Yscan","beam1RestVerticesUnfold_Yscan",xVar,yVar,x0_2,y0_2,w1N,w1M,rho_N1,xWidthN1,yWidthN1,rho_M1,xWidthM1,yWidthM1,rho_W1,xWidthW1,yWidthW1,w2N,w2M,xWidthN2,xWidthM2,xWidthW2,vtxRes);
  TripleGauss_V1  beam2RestVerticesUnfold_XScan("beam2RestVerticesUnfold_Xscan","beam2RestVerticesUnfold_Xscan",xVar,yVar,x0_12,y0_12,w2N,w2M,rho_N2,xWidthN2,yWidthN2,rho_M2,xWidthM2,yWidthM2,rho_W2,xWidthW2,yWidthW2,w1N,w1M,yWidthN1,xWidthM1,yWidthW1,vtxRes);
  TripleGauss_V2  beam2RestVerticesUnfold_YScan("beam2RestVerticesUnfold_Yscan","beam2RestVerticesUnfold_Yscan",xVar,yVar,x0_22,y0_22,w2N,w2M,rho_N2,xWidthN2,yWidthN2,rho_M2,xWidthM2,yWidthM2,rho_W2,xWidthW2,yWidthW2,w1N,w1M,xWidthN1,xWidthM1,xWidthW1,vtxRes);




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

  RooFitResult* r= simPdf.fitTo(combData,RooFit::PrintLevel(3),RooFit::Verbose(1),RooFit::Save());

  TF2 *multBeam = new TF2("multBeam",beamMultTG,-30,30,-30,30,26);
  multBeam->SetParameter(0, 0.0);
  multBeam->SetParameter(1, 0.0);
  multBeam->SetParameter(2,xWidthN1.getValV());
  multBeam->SetParameter(3,yWidthN1.getValV());
  multBeam->SetParameter(4,rho_N1.getValV());
  multBeam->SetParameter(5,xWidthM1.getValV());
  multBeam->SetParameter(6,yWidthM1.getValV());
  multBeam->SetParameter(7,rho_M1.getValV());
  multBeam->SetParameter(8,xWidthW1.getValV());
  multBeam->SetParameter(9,yWidthW1.getValV());
  multBeam->SetParameter(10,rho_W1.getValV());
  multBeam->SetParameter(11,w1N.getValV());
  multBeam->SetParameter(12,w1M.getValV());
  multBeam->SetParameter(13, 0.0);
  multBeam->SetParameter(14, 0.0);
  multBeam->SetParameter(15,xWidthN2.getValV());
  multBeam->SetParameter(16,yWidthN2.getValV());
  multBeam->SetParameter(17,rho_N2.getValV());
  multBeam->SetParameter(18,xWidthM2.getValV());
  multBeam->SetParameter(19,yWidthM2.getValV());
  multBeam->SetParameter(20,rho_M2.getValV());
  multBeam->SetParameter(21,xWidthW2.getValV());
  multBeam->SetParameter(22,yWidthW2.getValV());
  multBeam->SetParameter(23,rho_W2.getValV());
  multBeam->SetParameter(24,w2N.getValV());
  multBeam->SetParameter(25,w2M.getValV());

  std::cout<<"Overlap-Integral Fit: "<<multBeam->Integral(-30,30,-30,30)<<std::endl;

/*
   yWidth1N_fitf = yWidthN1.getValV();
   xWidth1N_fitf = xWidthN1.getValV();
   corr1N_fitf = rho_N1.getValV();
   yWidth1M_fitf = yWidthM1.getValV();
   xWidth1M_fitf = xWidthM1.getValV();
   corr1M_fitf = rho_M1.getValV();
   yWidth1W_fitf = yWidthW1.getValV();
   xWidth1W_fitf = xWidthW1.getValV();
   corr1W_fitf = rho_W1.getValV();
   yWidth2N_fitf = yWidthN2.getValV();
   xWidth2N_fitf = xWidthN2.getValV();
   corr2N_fitf = rho_N2.getValV();
   yWidth2M_fitf = yWidthM2.getValV();
   xWidth2M_fitf = xWidthM2.getValV();
   corr2M_fitf = rho_M2.getValV();
   yWidth2W_fitf = yWidthW2.getValV();
   xWidth2W_fitf = xWidthW2.getValV();
   corr2W_fitf = rho_W2.getValV();





float overlapReg = (reader->EvaluateRegression("MLP"))[0];
std::cout<<"Overlap-Integral Fit Regression: "<<overlapReg<<std::endl;
*/


    TRandom3 rand;
    rand.SetSeed(0);
    TFile *fAna = new TFile("DataAnalysisBunch"+bunchStr[i]+"TG_new_"+suff+".root","recreate");
    TH1F *integ = new TH1F("integ","Overlap Distribution",1000,0.,1.);

    TH1F *overlapInt_h = new TH1F("overlapInt","overlapInt",200,0.0,1.);
    overlapInt_h->Fill(multBeam->Integral(-30,30,-30,30));

    TH1F *xwidth1N_h = new TH1F("xwidth1N_h","xwidth1N_h",200,0.5,6.);
    xwidth1N_h->Fill(xWidthN1.getValV());
    TH1F *xwidth2N_h = new TH1F("xwidth2N_h","xwidth2N_h",200,0.5,6.);
    xwidth2N_h->Fill(xWidthN2.getValV());
    TH1F *ywidth1N_h = new TH1F("ywidth1N_h","ywidth1N_h",200,0.5,6.);
    ywidth1N_h->Fill(yWidthN1.getValV());
    TH1F *ywidth2N_h = new TH1F("ywidth2N_h","ywidth2N_h",200,0.5,6.);
    ywidth2N_h->Fill(yWidthN2.getValV());

    TH1F *xwidth1M_h = new TH1F("xwidth1M_h","xwidth1M_h",200,0.5,6.);
    xwidth1M_h->Fill(xWidthM1.getValV());
    TH1F *xwidth2M_h = new TH1F("xwidth2M_h","xwidth2M_h",200,0.5,6.);
    xwidth2M_h->Fill(xWidthM2.getValV());
    TH1F *ywidth1M_h = new TH1F("ywidth1M_h","ywidth1M_h",200,0.5,6.);
    ywidth1M_h->Fill(yWidthM1.getValV());
    TH1F *ywidth2M_h = new TH1F("ywidth2M_h","ywidth2M_h",200,0.5,6.);
    ywidth2M_h->Fill(yWidthM2.getValV());

    TH1F *xwidth1W_h = new TH1F("xwidth1W_h","xwidth1W_h",200,0.5,6.);
    xwidth1W_h->Fill(xWidthW1.getValV());
    TH1F *xwidth2W_h = new TH1F("xwidth2W_h","xwidth2W_h",200,0.5,6.);
    xwidth2W_h->Fill(xWidthW2.getValV());
    TH1F *ywidth1W_h = new TH1F("ywidth1W_h","ywidth1W_h",200,0.5,6.);
    ywidth1W_h->Fill(yWidthW1.getValV());
    TH1F *ywidth2W_h = new TH1F("ywidth2W_h","ywidth2W_h",200,0.5,6.);
    ywidth2W_h->Fill(yWidthW2.getValV());

    TH1F *weight1N_h = new TH1F("weight1N_h","weight1N_h",200,0.0,1.);
    weight1N_h->Fill(w1N.getValV());
    TH1F *weight2N_h = new TH1F("weight2N_h","weight2N_h",200,0.0,1.);
    weight2N_h->Fill(w2N.getValV());
    TH1F *weight1M_h = new TH1F("weight1M_h","weight1M_h",200,0.0,1.);
    weight1M_h->Fill(w1M.getValV());
    TH1F *weight2M_h = new TH1F("weight2M_h","weight2M_h",200,0.0,1.);
    weight2M_h->Fill(w2M.getValV());
    TH1F *weight1W_h = new TH1F("weight1W_h","weight1W_h",200,0.0,1.);
    weight1W_h->Fill(w1W.getValV());
    TH1F *weight2W_h = new TH1F("weight2W_h","weight2W_h",200,0.0,1.);
    weight2W_h->Fill(w2W.getValV());

    TH1F *rho1N_h = new TH1F("rho1N_h","rho1N_h",200,-0.5,0.5);
    rho1N_h->Fill(rho_N1.getValV());
    TH1F *rho2N_h = new TH1F("rho2N_h","rho2N_h",200,-0.5,0.5);
    rho2N_h->Fill(rho_N2.getValV());
    TH1F *rho1M_h = new TH1F("rho1M_h","rho1M_h",200,-0.5,0.5);
    rho1M_h->Fill(rho_M1.getValV());
    TH1F *rho2M_h = new TH1F("rho2M_h","rho2M_h",200,-0.5,0.5);
    rho2M_h->Fill(rho_M2.getValV());
    TH1F *rho1W_h = new TH1F("rho1W_h","rho1W_h",200,-0.5,0.5);
    rho1W_h->Fill(rho_W1.getValV());
    TH1F *rho2W_h = new TH1F("rho2W_h","rho2W_h",200,-0.5,0.5);
    rho2W_h->Fill(rho_W2.getValV());

    TH1F *y0_1_h = new TH1F("y0_1_h","y0_1_h",200,-0.5,0.5);
    y0_1_h->Fill(y0_1.getValV());
    TH1F *x0_1_h = new TH1F("x0_1_h","x0_1_h",200,-0.5,0.5);
    x0_1_h->Fill(x0_1.getValV());

    TH1F *y0_2_h = new TH1F("y0_2_h","y0_2_h",200,-0.5,0.5);
    y0_2_h->Fill(y0_2.getValV());
    TH1F *x0_2_h = new TH1F("x0_2_h","x0_2_h",200,-0.5,0.5);
    x0_2_h->Fill(x0_2.getValV());

    TH1F *y0_12_h = new TH1F("y0_12_h","y0_12_h",200,-0.5,0.5);
    y0_12_h->Fill(y0_12.getValV());
    TH1F *x0_12_h = new TH1F("x0_12_h","x0_12_h",200,-0.5,0.5);
    x0_12_h->Fill(x0_12.getValV());

    TH1F *y0_22_h = new TH1F("y0_22_h","y0_22_h",200,-0.5,0.5);
    y0_22_h->Fill(y0_22.getValV());
    TH1F *x0_22_h = new TH1F("x0_22_h","x0_22_h",200,-0.5,0.5);
    x0_22_h->Fill(x0_22.getValV());



 //Errors
    TH1F *xwidth1N_error_h = new TH1F("xwidth1N_error_h","xwidth1N_error_h",200,0.0,4.);
    xwidth1N_error_h->Fill(xWidthN1.getError());
    TH1F *xwidth2N_error_h = new TH1F("xwidth2N_error_h","xwidth2N_error_h",200,0.0,4.);
    xwidth2N_error_h->Fill(xWidthN2.getError());
    TH1F *ywidth1N_error_h = new TH1F("ywidth1N_error_h","ywidth1N_error_h",200,0.0,4.);
    ywidth1N_error_h->Fill(yWidthN1.getError());
    TH1F *ywidth2N_error_h = new TH1F("ywidth2N_error_h","ywidth2N_error_h",200,0.0,4.);
    ywidth2N_error_h->Fill(yWidthN2.getError());

    TH1F *xwidth1M_error_h = new TH1F("xwidth1M_error_h","xwidth1M_error_h",200,0.0,4.);
    //xwidth1M_error_h->Fill(xWidthM1.getError());
    double xwidth1M_error = TMath::Sqrt(TMath::Power(xWidthN1.getError(),2.0)+TMath::Power(xWidthM1Diff.getError(),2.0));
    xwidth1M_error_h->Fill(xwidth1M_error);
    TH1F *xwidth2M_error_h = new TH1F("xwidth2M_error_h","xwidth2M_error_h",200,0.0,4.);
    //xwidth2M_error_h->Fill(xWidthM2.getError());
    double xwidth2M_error = TMath::Sqrt(TMath::Power(xWidthN2.getError(),2.0)+TMath::Power(xWidthM2Diff.getError(),2.0));
    xwidth2M_error_h->Fill(xwidth2M_error);
    TH1F *ywidth1M_error_h = new TH1F("ywidth1M_error_h","ywidth1M_error_h",200,0.0,4.);
    //ywidth1M_error_h->Fill(yWidthM1.getError());
    double ywidth1M_error = TMath::Sqrt(TMath::Power(yWidthN1.getError(),2.0)+TMath::Power(xWidthM1Diff.getError(),2.0));
    ywidth1M_error_h->Fill(ywidth1M_error);
    TH1F *ywidth2M_error_h = new TH1F("ywidth2M_error_h","ywidth2M_error_h",200,0.0,4.);
    //ywidth2M_error_h->Fill(yWidthM2.getError());
    double ywidth2M_error = TMath::Sqrt(TMath::Power(yWidthN2.getError(),2.0)+TMath::Power(yWidthM2Diff.getError(),2.0));
    ywidth2M_error_h->Fill(ywidth2M_error);

    TH1F *xwidth1W_error_h = new TH1F("xwidth1W_error_h","xwidth1W_error_h",200,0.0,4.);
    //xwidth1W_error_h->Fill(xWidthW1.getError());
    double xwidth1W_error = TMath::Sqrt(TMath::Power(xWidthN1.getError(),2.0)+TMath::Power(xWidthM1Diff.getError(),2.0)+TMath::Power(xWidthW1Diff.getError(),2.0));
    xwidth1W_error_h->Fill(xwidth1W_error);
    TH1F *xwidth2W_error_h = new TH1F("xwidth2W_error_h","xwidth2W_error_h",200,0.0,4.);
    //xwidth2W_error_h->Fill(xWidthW2.getError());
    double xwidth2W_error = TMath::Sqrt(TMath::Power(xWidthN2.getError(),2.0)+TMath::Power(xWidthM2Diff.getError(),2.0)+TMath::Power(xWidthW2Diff.getError(),2.0));
    xwidth2W_error_h->Fill(xwidth2W_error);
    TH1F *ywidth1W_error_h = new TH1F("ywidth1W_error_h","ywidth1W_error_h",200,0.0,4.);
    //ywidth1W_error_h->Fill(yWidthW1.getError());
    double ywidth1W_error = TMath::Sqrt(TMath::Power(yWidthN1.getError(),2.0)+TMath::Power(yWidthM1Diff.getError(),2.0)+TMath::Power(yWidthW1Diff.getError(),2.0));
    ywidth1W_error_h->Fill(ywidth1W_error);
    TH1F *ywidth2W_error_h = new TH1F("ywidth2W_error_h","ywidth2W_error_h",200,0.0,4.);
    //ywidth2W_error_h->Fill(yWidthW2.getError());
    double ywidth2W_error = TMath::Sqrt(TMath::Power(yWidthN2.getError(),2.0)+TMath::Power(yWidthM2Diff.getError(),2.0)+TMath::Power(yWidthW2Diff.getError(),2.0));
    ywidth2W_error_h->Fill(ywidth2W_error);

    TH1F *weight1N_error_h = new TH1F("weight1N_error_h","weight1N_error_h",200,0.0,1.);
    double wN1_Error = 2.0 * sin(theta1.getValV()) * cos(phi1.getValV()) * TMath::Sqrt(TMath::Power(theta1.getError()*cos(theta1.getValV())*cos(phi1.getValV()), 2.0) + TMath::Power(phi1.getError()*sin(theta1.getValV())*sin(phi1.getValV()), 2.0));
    weight1N_error_h->Fill(wN1_Error);
    TH1F *weight1M_error_h = new TH1F("weight1M_error_h","weight1M_error_h",200,0.0,1.);
    double wM1_Error = 2.0 * sin(theta1.getValV()) * sin(phi1.getValV()) * TMath::Sqrt(TMath::Power(theta1.getError()*cos(theta1.getValV())*sin(phi1.getValV()), 2.0) + TMath::Power(phi1.getError()*sin(theta1.getValV())*cos(phi1.getValV()), 2.0));
    weight1M_error_h->Fill(wM1_Error);
    TH1F *weight1W_error_h = new TH1F("weight1W_error_h","weight1W_error_h",200,0.0,1.);
    double wW1_Error = 2.0 * sin(theta1.getValV()) * cos(theta1.getValV()) * theta1.getError();
    weight1W_error_h->Fill(wW1_Error);
    TH1F *weight2N_error_h = new TH1F("weight2N_error_h","weight2N_error_h",200,0.0,1.);
    double wN2_Error = 2.0 * sin(theta2.getValV()) * cos(phi2.getValV()) * TMath::Sqrt(TMath::Power(theta2.getError()*cos(theta2.getValV())*cos(phi2.getValV()), 2.0) + TMath::Power(phi2.getError()*sin(theta2.getValV())*sin(phi2.getValV()), 2.0));
    weight2N_error_h->Fill(wN2_Error);
    TH1F *weight2M_error_h = new TH1F("weight2M_error_h","weight2M_error_h",200,0.0,1.);
    double wM2_Error = 2.0 * sin(theta2.getValV()) * sin(phi2.getValV()) * TMath::Sqrt(TMath::Power(theta2.getError()*cos(theta2.getValV())*sin(phi2.getValV()), 2.0) + TMath::Power(phi2.getError()*sin(theta2.getValV())*cos(phi2.getValV()), 2.0));
    weight2M_error_h->Fill(wM2_Error);
    TH1F *weight2W_error_h = new TH1F("weight2W_error_h","weight2W_error_h",200,0.0,1.);
    double wW2_Error = 2.0 * sin(theta2.getValV()) * cos(theta2.getValV()) * theta2.getError();
    weight2W_error_h->Fill(wW2_Error);

    TH1F *rho1N_error_h = new TH1F("rho1N_error_h","rho1N_error_h",200,-0.5,0.5);
    rho1N_error_h->Fill(rho_N1.getError());
    TH1F *rho2N_error_h = new TH1F("rho2N_error_h","rho2N_error_h",200,-0.5,0.5);
    rho2N_error_h->Fill(rho_N2.getError());
    TH1F *rho1M_error_h = new TH1F("rho1M_error_h","rho1M_error_h",200,-0.5,0.5);
    rho1M_error_h->Fill(rho_M1.getError());
    TH1F *rho2M_error_h = new TH1F("rho2M_error_h","rho2M_error_h",200,-0.5,0.5);
    rho2M_error_h->Fill(rho_M2.getError());
    TH1F *rho1W_error_h = new TH1F("rho1W_error_h","rho1W_error_h",200,-0.5,0.5);
    rho1W_error_h->Fill(rho_W1.getError());
    TH1F *rho2W_error_h = new TH1F("rho2W_error_h","rho2W_error_h",200,-0.5,0.5);
    rho2W_error_h->Fill(rho_W2.getError());


  for(int k=0; k<1000;k++)
    {
      //rand.Gaus(vtxX,vertexResolution);
      multBeam->SetParameter(0, 0.0);
      multBeam->SetParameter(1, 0.0);
      multBeam->SetParameter(2,rand.Gaus(xWidthN1.getValV(),xWidthN1.getError()));
      multBeam->SetParameter(3,rand.Gaus(yWidthN1.getValV(),yWidthN1.getError()));
      multBeam->SetParameter(4,rand.Gaus(rho_N1.getValV(),rho_N1.getError()));
      //multBeam->SetParameter(5,rand.Gaus(xWidthM1.getValV(),xWidthM1.getError()));
      //multBeam->SetParameter(6,rand.Gaus(yWidthM1.getValV(),yWidthM1.getError()));
      multBeam->SetParameter(5,rand.Gaus(xWidthM1.getValV(),xwidth1M_error));
      multBeam->SetParameter(6,rand.Gaus(yWidthM1.getValV(),ywidth1M_error));
      multBeam->SetParameter(7,rand.Gaus(rho_M1.getValV(),rho_M1.getError()));
      //multBeam->SetParameter(8,rand.Gaus(xWidthW1.getValV(),xWidthW1.getError()));
      //multBeam->SetParameter(9,rand.Gaus(yWidthW1.getValV(),yWidthW1.getError()));
      multBeam->SetParameter(8,rand.Gaus(xWidthW1.getValV(),xwidth1W_error));
      multBeam->SetParameter(9,rand.Gaus(yWidthW1.getValV(),ywidth1W_error));
      multBeam->SetParameter(10,rand.Gaus(rho_W1.getValV(),rho_W1.getError()));
      double newtheta1 = rand.Gaus(theta1.getValV(),theta1.getError());
      double newphi1 = rand.Gaus(phi1.getValV(),phi1.getError());
      multBeam->SetParameter(11,TMath::Power(sin(newtheta1)*cos(newphi1),2.0));
      multBeam->SetParameter(12,TMath::Power(sin(newtheta1)*sin(newphi1),2.0));
      multBeam->SetParameter(13, 0.0);
      multBeam->SetParameter(14, 0.0);
      multBeam->SetParameter(15,rand.Gaus(xWidthN2.getValV(),xWidthN2.getError()));
      multBeam->SetParameter(16,rand.Gaus(yWidthN2.getValV(),yWidthN2.getError()));
      multBeam->SetParameter(17,rand.Gaus(rho_N2.getValV(),rho_N2.getError()));
      //multBeam->SetParameter(18,rand.Gaus(xWidthM2.getValV(),xWidthM2.getError()));
      //multBeam->SetParameter(19,rand.Gaus(yWidthM2.getValV(),yWidthM2.getError()));
      multBeam->SetParameter(18,rand.Gaus(xWidthM2.getValV(),xwidth2M_error));
      multBeam->SetParameter(19,rand.Gaus(yWidthM2.getValV(),ywidth2M_error));
      multBeam->SetParameter(20,rand.Gaus(rho_M2.getValV(),rho_M2.getError()));
      //multBeam->SetParameter(21,rand.Gaus(xWidthW2.getValV(),xWidthW2.getError()));
      //multBeam->SetParameter(22,rand.Gaus(yWidthW2.getValV(),yWidthW2.getError()));
      multBeam->SetParameter(21,rand.Gaus(xWidthW2.getValV(),xwidth2W_error));
      multBeam->SetParameter(22,rand.Gaus(yWidthW2.getValV(),ywidth2W_error));
      multBeam->SetParameter(23,rand.Gaus(rho_W2.getValV(),rho_W2.getError()));
      double newtheta2 = rand.Gaus(theta2.getValV(),theta2.getError());
      double newphi2 = rand.Gaus(phi2.getValV(),phi2.getError());
      multBeam->SetParameter(24,TMath::Power(sin(newtheta2)*cos(newphi2),2.0));
      multBeam->SetParameter(25,TMath::Power(sin(newtheta2)*sin(newphi2),2.0));

      cout<<multBeam->Integral(-30,30,-30,30)<<endl;

      integ->Fill(multBeam->Integral(-30,30,-30,30));
    }



  int nbins = 5*19;

 TH2D* hmodelX1 = (TH2D*) beam1RestVerticesUnfold_XScan.createHistogram("hmodelX1",xVar,Binning(nbins),YVar(yVar,Binning(nbins)))  ;
 //RooDataHist* hmodelX1_h =  beam1RestVerticesUnfold_XScan.generateBinned(RooArgSet(xVar,yVar),1000000);
 //TH2D* hmodelX1 = (TH2D*) hmodelX1_h->createHistogram("hmodelX1",xVar,Binning(19),YVar(yVar,Binning(19)))  ;
 TH2D* hmodelY1 = (TH2D*) beam1RestVerticesUnfold_YScan.createHistogram("hmodelY1",xVar,Binning(nbins),YVar(yVar,Binning(nbins)))  ;
 TH2D* hmodelX2 = (TH2D*) beam2RestVerticesUnfold_XScan.createHistogram("hmodelX2",xVar,Binning(nbins),YVar(yVar,Binning(nbins)))  ;
 TH2D* hmodelY2 = (TH2D*) beam2RestVerticesUnfold_YScan.createHistogram("hmodelY2",xVar,Binning(nbins),YVar(yVar,Binning(nbins)))  ;

  TH2D* hdataX1 = (TH2D*) scanXBeam1RestDataHist->createHistogram("hdataX1",xVar,Binning(nbins),YVar(yVar,Binning(nbins))) ;
  TH2D* hdataY1 = (TH2D*) scanYBeam1RestDataHist->createHistogram("hdataY1",xVar,Binning(nbins),YVar(yVar,Binning(nbins))) ;
  TH2D* hdataX2 = (TH2D*) scanXBeam2RestDataHist->createHistogram("hdataX2",xVar,Binning(nbins),YVar(yVar,Binning(nbins))) ;
  TH2D* hdataY2 = (TH2D*) scanYBeam2RestDataHist->createHistogram("hdataY2",xVar,Binning(nbins),YVar(yVar,Binning(nbins))) ;

  hmodelX1->Scale(hdataX1->Integral());
  hmodelY1->Scale(hdataY1->Integral());
  hmodelX2->Scale(hdataX2->Integral());
  hmodelY2->Scale(hdataY2->Integral());


 integ->Write("Fit_Stat_error");

  /*TH2D* resX1 = (TH2D*) hdataX1->Clone();
  TH2D* resY1 = (TH2D*) hdataY1->Clone();
  TH2D* resX2 = (TH2D*) hdataX2->Clone();
  TH2D* resY2 = (TH2D*) hdataY2->Clone();*/

  TH2D* resX1 = new TH2D("BeamImageX1","BeamImageX1",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* resY1 = new TH2D("BeamImageY1","BeamImageY1",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* resX2 = new TH2D("BeamImageX2","BeamImageX2",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* resY2 = new TH2D("BeamImageY2","BeamImageY2",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);

  TH2D* dataHistX1 = new TH2D("dataHistX1","dataHistX1",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* dataHistY1 = new TH2D("dataHistY1","dataHistY1",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* dataHistX2 = new TH2D("dataHistX2","dataHistX2",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* dataHistY2 = new TH2D("dataHistY2","dataHistY2",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);

  TH2D* modelHistX1 = new TH2D("modelHistX1","modelHistX1",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* modelHistY1 = new TH2D("modelHistY1","modelHistY1",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* modelHistX2 = new TH2D("modelHistX2","modelHistX2",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);
  TH2D* modelHistY2 = new TH2D("modelHistY2","modelHistY2",nbins,-10*scaling,10*scaling,nbins,-10*scaling,10*scaling);

  double chi2Y1 = 0.;
  double chi2Y2 = 0.;
  double chi2X1 = 0.;
  double chi2X2 = 0.;

  double dofX1 = 0.;
  double dofY1 = 0.;
  double dofX2 = 0.;
  double dofY2 = 0.;


  for(int xs=0; xs<hdataX1->GetXaxis()->GetNbins(); xs++)
    {
      for(int ys=0; ys<hdataX1->GetYaxis()->GetNbins(); ys++)
	{
	  double valDiffX1;
	  double statErrX1;
	  if( hdataX1->GetBinError(xs,ys)>0. &&  hdataX1->GetBinContent(xs,ys)>0.){
	    dataHistX1->SetBinContent(xs,ys,hdataX1->GetBinContent(xs,ys));
	    modelHistX1->SetBinContent(xs,ys,hmodelX1->GetBinContent(xs,ys));
	  valDiffX1 = hdataX1->GetBinContent(xs,ys) - hmodelX1->GetBinContent(xs,ys);
	  statErrX1 = hdataX1->GetBinError(xs,ys);
	  resX1->SetBinContent(xs,ys,valDiffX1/statErrX1);
	  dofX1 +=1.;
	  chi2X1 +=(TMath::Power(valDiffX1/statErrX1,2.));///hmodelX1->GetBinContent(xs,ys);

	  }
	  //else{resX1->SetBinContent(xs,ys,0.);}
	  double valDiffY1;
	  double statErrY1;
	  if( hdataY1->GetBinError(xs,ys)>0. &&  hdataY1->GetBinContent(xs,ys)>0.){
	    dataHistY1->SetBinContent(xs,ys,hdataY1->GetBinContent(xs,ys));
	    modelHistY1->SetBinContent(xs,ys,hmodelY1->GetBinContent(xs,ys));
	  valDiffY1 = hdataY1->GetBinContent(xs,ys) - hmodelY1->GetBinContent(xs,ys);
	  statErrY1 = hdataY1->GetBinError(xs,ys);
	  resY1->SetBinContent(xs,ys,valDiffY1/statErrY1);
	  chi2Y1+=(TMath::Power(valDiffY1/statErrY1,2.));
	  dofY1 +=1.;
	  }
	  //else{resY1->SetBinContent(xs,ys,0.);}
	  double valDiffX2;
	  double statErrX2;
	  if( hdataX2->GetBinError(xs,ys)>0. &&  hdataX2->GetBinContent(xs,ys)>0.){
	    dataHistX2->SetBinContent(xs,ys,hdataX2->GetBinContent(xs,ys));
	    modelHistX2->SetBinContent(xs,ys,hmodelX2->GetBinContent(xs,ys));
	  valDiffX2 = hdataX2->GetBinContent(xs,ys) - hmodelX2->GetBinContent(xs,ys);
	  statErrX2 = hdataX2->GetBinError(xs,ys);
	  resX2->SetBinContent(xs,ys,valDiffX2/statErrX2);
	  chi2X2+=(TMath::Power(valDiffX2/statErrX2,2.));
	  dofX2 +=1.;
	  }
	  else{resX2->SetBinContent(xs,ys,0.);}
	  double valDiffY2;
	  double statErrY2;
	  if( hdataY2->GetBinError(xs,ys)>0. &&  hdataY2->GetBinContent(xs,ys)>0.){
	    dataHistY2->SetBinContent(xs,ys,hdataY2->GetBinContent(xs,ys));
	    modelHistY2->SetBinContent(xs,ys,hmodelY2->GetBinContent(xs,ys));
	  valDiffY2 = hdataY2->GetBinContent(xs,ys) - hmodelY2->GetBinContent(xs,ys);
	  statErrY2 = hdataY2->GetBinError(xs,ys);
	  resY2->SetBinContent(xs,ys,valDiffY2/statErrY2);
	  chi2Y2+=(TMath::Power(valDiffY2/statErrY2,2.));
	  dofY2 +=1.;
	  }
	  else{resY2->SetBinContent(xs,ys,0.);}
	}
    }

  // cout<<hdata->GetYaxis()->GetNbins()<<endl;
  // cout<<hdata->GetBinError(10,10)<<endl;

  cout<<"chi2 X1: "<<(chi2X1)<<endl;
  cout<<"dof X1: "<<dofX1<<endl;
cout<<"chi2 X2: "<<(chi2X2)<<endl;
  cout<<"dof X2: "<<dofX2<<endl;
cout<<"chi2 Y1: "<<(chi2Y1)<<endl;
  cout<<"dof Y1: "<<dofY1<<endl;
cout<<"chi2 Y2: "<<(chi2Y2)<<endl;
  cout<<"dof Y2: "<<dofY2<<endl;


  TCanvas *cz = new TCanvas("cz","cz",400,600);
  resX1->SetOption("Colz");
  resY1->SetOption("Colz");
  resX2->SetOption("Colz");
  resY2->SetOption("Colz");
  resX1->SetTitle("Beam Image X1");
  resY1->SetTitle("Beam Image Y1");
  resX2->SetTitle("Beam Image X2");
  resY2->SetTitle("Beam Image Y2");

  resX1->GetXaxis()->SetTitle("[cm]");
  resX1->GetYaxis()->SetTitle("[cm]");
  resX2->GetXaxis()->SetTitle("[cm]");
  resX2->GetYaxis()->SetTitle("[cm]");

  resY1->GetXaxis()->SetTitle("[cm]");
  resY1->GetYaxis()->SetTitle("[cm]");
  resY2->GetXaxis()->SetTitle("[cm]");
  resY2->GetYaxis()->SetTitle("[cm]");


  resX1->Write("resX1");
  resY1->Write("resY1");
  resX2->Write("resX2");
  resY2->Write("resY2");

  dataHistX1->Write("dataHistX1");
  dataHistY1->Write("dataHistY1");
  dataHistX2->Write("dataHistX2");
  dataHistY2->Write("dataHistY2");

  modelHistX1->Write("modelHistX1");
  modelHistY1->Write("modelHistY1");
  modelHistX2->Write("modelHistX2");
  modelHistY2->Write("modelHistY2");

  r->Write("fitResult");

  fAna->Write();

 //}
}
