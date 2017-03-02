//#include "MyPdfV1.h"
#include "../../beamdensities/SingleGauss_V1.h"
#include "../../beamdensities/SingleGauss_V2.h"
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


Double_t beamMultSGt(Double_t *x,Double_t *par)
{
  Double_t arg = 0;
  Double_t pi = 3.1415926;

  double x01 = par[0];
  double y01 = par[1];
  double xwidthN1 = par[2];
  double ywidthN1 = par[3];
  double rhoN1 = par[4];

  double x02 = par[5];
  double y02 = par[6];
  double xwidthN2 = par[7];
  double ywidthN2 = par[8];
  double rhoN2 = par[9];

  double beamN1_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN1,2)))*TMath::Sqrt(TMath::Power(xwidthN1,2.))*TMath::Sqrt(TMath::Power(ywidthN1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN1,2))*(TMath::Power((x[0]-x01)/xwidthN1,2.0)+TMath::Power((x[1]-y01)/ywidthN1,2.0)-2*rhoN1*(x[0]-x01)*(x[1]-y01)/(xwidthN1*ywidthN1)));
  double beamN2_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN2,2)))*TMath::Sqrt(TMath::Power(xwidthN2,2.))*TMath::Sqrt(TMath::Power(ywidthN2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN2,2))*(TMath::Power((x[0]-x02)/xwidthN2,2.0)+TMath::Power((x[1]-y02)/ywidthN2,2.0)-2*rhoN2*(x[0]-x02)*(x[1]-y02)/(xwidthN2*ywidthN2)));

  double product2 = beamN1_ * beamN2_;

  return product2;


}


void vdmScanTreeAnalyzerSG_onData2_studyUnc_2016()
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
TString bunchStr[5] = {"41","281","872","1783","2063"};
 gStyle->SetOptStat(0);


/*TMVA::Tools::Instance();

TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );



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


 for(int i=0;i<5;i++){

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

  RooRealVar yWidthN1("yWidthN1","yWidthN1",1.3,3.0) ;
  RooRealVar xWidthN1("xWidthN1","xWidthN1",1.3,3.0) ;
  RooRealVar rho_N1("rho_N1","rho_N1",-0.48,0.48) ;

  RooRealVar yWidthN2("yWidthN2","yWidthN2",1.3,3.0) ;
  RooRealVar xWidthN2("xWidthN2","xWidthN2",1.3,3.0) ;
  RooRealVar rho_N2("rho_N2","rho_N2",-0.48,0.48) ;

  RooRealVar vtxRes("vtxRes","vtxRes",0.00356/scaling) ;
  vtxRes.setConstant();

  RooRealVar xWidthRes("xWidthRes","xWidthRes",0.62,0.62) ;
  RooRealVar yWidthRes("yWidthRes","yWidthRes",0.62,0.62) ;
  RooRealVar meanRes("meanRes","meanRes",0.0,0.0) ;



  RooGaussian resX("resX","resX",xVar,meanRes,xWidthRes);
  RooGaussian resY("resY","resY",yVar,meanRes,yWidthRes);


  SingleGauss_V1  beam1RestVerticesUnfold_XScan("beam1RestVerticesUnfold_Xscan","beam1RestVerticesUnfold_Xscan",xVar,yVar,x0_1,y0_1,rho_N1,xWidthN1,yWidthN1,yWidthN2,vtxRes);
  SingleGauss_V2  beam1RestVerticesUnfold_YScan("beam1RestVerticesUnfold_Yscan","beam1RestVerticesUnfold_Yscan",xVar,yVar,x0_2,y0_2,rho_N1,xWidthN1,yWidthN1,xWidthN2,vtxRes);
  SingleGauss_V1  beam2RestVerticesUnfold_XScan("beam2RestVerticesUnfold_Xscan","beam2RestVerticesUnfold_Xscan",xVar,yVar,x0_12,y0_12,rho_N2,xWidthN2,yWidthN2,yWidthN1,vtxRes);
  SingleGauss_V2  beam2RestVerticesUnfold_YScan("beam2RestVerticesUnfold_Yscan","beam2RestVerticesUnfold_Yscan",xVar,yVar,x0_22,y0_22,rho_N2,xWidthN2,yWidthN2,xWidthN1,vtxRes);




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



TF2 *multBeam = new TF2("multBeam",beamMultSGt,-30,30,-30,30,10);




  multBeam->SetParameter(0, 0.0);
  multBeam->SetParameter(1, 0.0);
  multBeam->SetParameter(2,xWidthN1.getValV());
  multBeam->SetParameter(3,yWidthN1.getValV());
  multBeam->SetParameter(4,rho_N1.getValV());
  multBeam->SetParameter(5, 0.0);
  multBeam->SetParameter(6, 0.0);
  multBeam->SetParameter(7,xWidthN2.getValV());
  multBeam->SetParameter(8,yWidthN2.getValV());
  multBeam->SetParameter(9,rho_N2.getValV());

  std::cout<<"Overlap-Integral Fit: "<<multBeam->Integral(-30,30,-30,30)<<std::endl;

/*
  yWidth1N_fitf = yWidthN1.getValV();
  xWidth1N_fitf = xWidthN1.getValV();
  corr1N_fitf = rho_N1.getValV();
  yWidth1W_fitf = yWidthW1.getValV();
  xWidth1W_fitf = xWidthW1.getValV();
  corr1W_fitf = rho_W1.getValV();
   weight1_fitf = w1.getValV();
   yWidth2N_fitf = yWidthN2.getValV();
   xWidth2N_fitf = xWidthN2.getValV();
   corr2N_fitf = rho_N2.getValV();
   yWidth2W_fitf = yWidthW2.getValV();
   xWidth2W_fitf = xWidthW2.getValV();
   corr2W_fitf = rho_W2.getValV();
   weight2_fitf = w2.getValV();






float overlapReg = (reader->EvaluateRegression("MLP"))[0];

std::cout<<"Overlap-Integral Fit Regression: "<<overlapReg<<std::endl;
*/


TRandom3 rand;
  rand.SetSeed(0);
 TFile *fAna = new TFile("DataAnalysisBunch"+bunchStr[i]+"SG_new_"+suff+".root","recreate");
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

   TH1F *rho1N_h = new TH1F("rho1N_h","rho1N_h",200,-0.5,0.5);
   rho1N_h->Fill(rho_N1.getValV());
   TH1F *rho2N_h = new TH1F("rho2N_h","rho2N_h",200,-0.5,0.5);
    rho2N_h->Fill(rho_N2.getValV());

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

   TH1F *rho1N_error_h = new TH1F("rho1N_error_h","rho1N_error_h",200,-0.5,0.5);
   rho1N_error_h->Fill(rho_N1.getError());
   TH1F *rho2N_error_h = new TH1F("rho2N_error_h","rho2N_error_h",200,-0.5,0.5);
    rho2N_error_h->Fill(rho_N2.getError());

  for(int k=0; k<1000;k++)
    {
      //rand.Gaus(vtxX,vertexResolution);
      multBeam->SetParameter(0, 0.0);
      multBeam->SetParameter(1, 0.0);
      multBeam->SetParameter(2,rand.Gaus(xWidthN1.getValV(),xWidthN1.getError()));
      multBeam->SetParameter(3,rand.Gaus(yWidthN1.getValV(),yWidthN1.getError()));
      multBeam->SetParameter(4,rand.Gaus(rho_N1.getValV(),rho_N1.getError()));
      multBeam->SetParameter(5, 0.0);
      multBeam->SetParameter(6, 0.0);
      multBeam->SetParameter(7,rand.Gaus(xWidthN2.getValV(),xWidthN2.getError()));
      multBeam->SetParameter(8,rand.Gaus(yWidthN2.getValV(),yWidthN2.getError()));
      multBeam->SetParameter(9,rand.Gaus(rho_N2.getValV(),rho_N2.getError()));

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

 }
}
