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
//#include "MyPdfV3.h"
//#include "MyPdfV4.h"
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
#include "RooHist.h"
#include "TStyle.h"
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
  
  double product = (beamN1+beamW1)*(beamN2+beamW2);

  double beamN1_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN1,2)))*TMath::Sqrt(TMath::Power(xwidthN1,2.))*TMath::Sqrt(TMath::Power(ywidthN1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN1,2))*(TMath::Power((x[0]-x01)/xwidthN1,2.0)+TMath::Power((x[1]-y01)/ywidthN1,2.0)-2*rhoN1*(x[0]-x01)*(x[1]-y01)/(xwidthN1*ywidthN1)));

  double beamW1_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW1,2)))*TMath::Sqrt(TMath::Power(xwidthW1,2.))*TMath::Sqrt(TMath::Power(ywidthW1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW1,2))*(TMath::Power((x[0]-x01)/xwidthW1,2.0)+TMath::Power((x[1]-y01)/ywidthW1,2.0)-2*rhoW1*(x[0]-x01)*(x[1]-y01)/(xwidthW1*ywidthW1)));

  double beamN2_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN2,2)))*TMath::Sqrt(TMath::Power(xwidthN2,2.))*TMath::Sqrt(TMath::Power(ywidthN2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN2,2))*(TMath::Power((x[0]-x02)/xwidthN2,2.0)+TMath::Power((x[1]-y02)/ywidthN2,2.0)-2*rhoN2*(x[0]-x02)*(x[1]-y02)/(xwidthN2*ywidthN2)));

  double beamW2_ = 1./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW2,2)))*TMath::Sqrt(TMath::Power(xwidthW2,2.))*TMath::Sqrt(TMath::Power(ywidthW2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW2,2))*(TMath::Power((x[0]-x02)/xwidthW2,2.0)+TMath::Power((x[1]-y02)/ywidthW2,2.0)-2*rhoW2*(x[0]-x02)*(x[1]-y02)/(xwidthW2*ywidthW2)));

  double product2 = (nw_weight1*beamN1_ + (1.-nw_weight1)*beamW1_) * (nw_weight2*beamN2_ + (1.-nw_weight2)*beamW2_);

  return product2;

}

void drawFit(RooRealVar vtxframevar,RooPlot* vtxframe,TCanvas* c,const char* name) {
  // Pulls
  RooHist* hpull = vtxframe->pullHist();
  RooPlot* pullframe = vtxframevar.frame();
  pullframe->SetTitle("  ");
  pullframe->GetXaxis()->SetTitleSize(0.09);
  pullframe->GetYaxis()->SetRangeUser(-4.5,4.5);
  pullframe->addPlotable(hpull,"P");
  c->Divide(1,2);
  c->cd(2)->SetPad(0,0,1.0,0.3); c->cd(2)->SetTopMargin(0.05); c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15); c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetGridy(1);
  pullframe->Draw();
  c->cd(1)->SetPad(0,0.3,1.0,1.0); c->cd(1)->SetTopMargin(0.1); c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15); c->cd(1)->SetRightMargin(0.07);
  vtxframe->Draw();
  c->Print(name);
  c->Close();
}

double findMax_v0(RooCurve* fitCurve,double mu1,double mu2) {
  // find maximum value of DG fit
  int delta = TMath::CeilNint((mu1-mu2)/100);
  double maxfitVal = 0;
  for(int m=0;m<100;m++) {
    double fitVal = fitCurve->Eval(TMath::Min(mu1,mu2)+m*delta);
    if(fitVal>maxfitVal) maxfitVal = fitVal;
  }
  return maxfitVal;
}

double findMax(RooCurve* fitCurve,double x1,double x2) {
  // find maximum value of DG fit
  double delta = 0.001;//TMath::CeilNint((x1-x2)/100);
  int nsteps = TMath::FloorNint(fabs(x1-x2)/0.001);
  double maxfitVal = 0;
  for(int m=0;m<nsteps;m++) {
    double fitVal = fitCurve->Eval(TMath::Min(x1,x2)+m*delta);
    if(fitVal>maxfitVal) maxfitVal = fitVal;
  }
  return maxfitVal;
}

void overlapDiff_DG(TString suffix)
{
  TRandom3 r;
  gStyle->SetOptStat(0);
  r.SetSeed(0);
  int scansteps = 25.;
TFile *fAna = new TFile("overlapDiff_DG_TOYS_2016"+suffix+".root","recreate");
   TTree *outTrainTree = new TTree("outTrainTree","tree of VdM fit parameters");

   TH1F *overDiff = new TH1F("overDiff","VdM Scan Measured vs True Overlap",200,-0.1,0.1);

 TH1F *capSigAeffDiff = new TH1F("capSigAeffDiff","VdM Scan Measured vs True Overlap",360,0.,0.12);

double in_yWidth1N ;
 double in_xWidth1N ;
double in_corr1N ;
double in_yWidth1W ;
 double in_xWidth1W ;
 double in_corr1W ;
 double in_weight1 ;
 double in_yWidth2N ;
 double in_xWidth2N;
 double in_corr2N ;
 double in_yWidth2W ;
 double in_xWidth2W ;
 double in_corr2W ;
 double in_weight2 ;
double in_overlapTrueTree;

//TFile *f = TFile::Open("DataAnalysisBunch"+suffix+"DG_with_Improved_addOrigin.root");

//TFile *f = TFile::Open("DataAnalysisBunch"+suffix+"DG_with_Improved.root");

//TFile *f = TFile::Open("DataAnalysisBunch"+suffix+"DG_with_Improved_addOriginPromptReco2015.root");

//TFile *f = TFile::Open("DataAnalysisBunch"+suffix+"DG_with_Improved_addOrigin.root");
//TFile *f = TFile::Open("DataAnalysisBunch2211DG_with_Improved_addOriginPromptReco2015.root");
TFile *f = TFile::Open("DataAnalysisBunch"+suffix+"DG_withML_StronRescale_NewScale.root");

   TFile *file = TFile::Open("../ultimate.root");
   TTree *tree = new TTree("outTrainTree", "outTrainTree");
   tree->SetBranchAddress("yWidth1N_true",&in_yWidth1N );
   tree->SetBranchAddress("xWidth1N_true",&in_xWidth1N );
   tree->SetBranchAddress("corr1N_true",&in_corr1N );
   tree->SetBranchAddress("yWidth1W_true",&in_yWidth1W );
   tree->SetBranchAddress("xWidth1W_true",&in_xWidth1W );
   tree->SetBranchAddress("corr1W_true",&in_corr1W );
   tree->SetBranchAddress("weight1_true",&in_weight1 );

   tree->SetBranchAddress("yWidth2N_true",&in_yWidth2N );
   tree->SetBranchAddress("xWidth2N_true",&in_xWidth2N );
   tree->SetBranchAddress("corr2N_true",&in_corr2N );
   tree->SetBranchAddress("yWidth2W_true",&in_yWidth2W );
   tree->SetBranchAddress("xWidth2W_true",&in_xWidth2W );
   tree->SetBranchAddress("corr2W_true",&in_corr2W );
   tree->SetBranchAddress("weight2_true",&in_weight2 );
tree->SetBranchAddress("overlapFit",&in_overlapTrueTree );
   for(int z = 0 ; z < 1 ; z++){

     // tree->GetEntry(z);
     //cout<<"Iteration "<<z<<endl;

   
  //51
  /*double in_yWidth1N = 1.77;//1.8+r.Uniform(0.7);
  double in_xWidth1N = 1.98;//1.8+r.Uniform(0.7);
  double in_corr1N = -0.101;//-0.4+r.Uniform(0.8);
  double in_yWidth1W = 2.29;//1.8+r.Uniform(0.7);
  double in_xWidth1W = 2.4;//1.8+r.Uniform(0.7);
  double in_corr1W = -0.33;//-0.4+r.Uniform(0.8);
  double in_weight1 = 0.785;//r.Uniform(1.);

  double in_yWidth2N = 1.74;//1.8+r.Uniform(0.7);
  double in_xWidth2N = 1.93;//1.8+r.Uniform(0.7);
  double in_corr2N = -0.184;//-0.4+r.Uniform(0.8);
    // cout<<in_corr2N<<endl;
    // cout<<in_corr1N<<endl;
  double in_yWidth2W = 2.3;//1.8+r.Uniform(0.7);
  double in_xWidth2W = 2.29;//1.8+r.Uniform(0.7);
  double in_corr2W = 0.3;//-0.4+r.Uniform(0.8);
  double in_weight2 = 0.712;//r.Uniform(1.);*/



  /*double in_yWidth1N = 1.719;//1.8+r.Uniform(0.7);
  double in_xWidth1N = 1.987;//1.8+r.Uniform(0.7);
  double in_corr1N = -0.0813;//-0.4+r.Uniform(0.8);
  double in_yWidth1W = 2.18;//1.8+r.Uniform(0.7);
  double in_xWidth1W = 2.26;//1.8+r.Uniform(0.7);
  double in_corr1W = -0.305;//-0.4+r.Uniform(0.8);
  double in_weight1 = 0.71;//r.Uniform(1.);

  double in_yWidth2N = 1.733;//1.8+r.Uniform(0.7);
  double in_xWidth2N = 1.93;//1.8+r.Uniform(0.7);
  double in_corr2N = -0.206;//-0.4+r.Uniform(0.8);
  double in_yWidth2W = 2.216;//1.8+r.Uniform(0.7);
  double in_xWidth2W = 2.23;//1.8+r.Uniform(0.7);
  double in_corr2W = 0.27;//-0.4+r.Uniform(0.8);
  double in_weight2 = 0.675;//r.Uniform(1.);
*/
  /*
   double in_yWidth1N = 1.6+r.Uniform(0.2);
  double in_xWidth1N = 1.6+r.Uniform(0.2);
  double in_corr1N = 0.2;//-0.4+r.Uniform(0.8);
  double in_yWidth1W = 2.4+r.Uniform(0.2);
  double in_xWidth1W = 2.4+r.Uniform(0.2);
  double in_corr1W = 0.2;//-0.4+r.Uniform(0.8);
  double in_weight1 = 0.5;//r.Uniform(1.);

  double in_yWidth2N = 1.6+r.Uniform(0.2);
  double in_xWidth2N = 1.6+r.Uniform(0.2);
  double in_corr2N = 0.2;//-0.4+r.Uniform(0.8);
  double in_yWidth2W =2.4+r.Uniform(0.2);
  double in_xWidth2W = 2.4+r.Uniform(0.2);
  double in_corr2W = 0.2;//-0.4+r.Uniform(0.8);
  double in_weight2 = 0.5;//r.Uniform(1.);*/

  /*double in_yWidth1N = 1.95464;
double in_xWidth1N = 1.87367;
double in_corr1N = 0.394601;
double in_yWidth1W = 2.01402;
double in_xWidth1W = 2.3581;
double in_corr1W = 0.120333;
double in_weight1 = 0.521179;
double in_yWidth2N = 1.93209;
double in_xWidth2N= 1.91124;
double in_corr2N = 0.0631183;
double in_yWidth2W = 2.37724;
double in_xWidth2W = 2.23537;
double in_corr2W = 0.302163;
double in_weight2 = 0.874279;*/

  /*double in_yWidth1N = 1.981;
double in_xWidth1N = 1.881;
double in_corr1N = 0.419 ;
double in_yWidth1W = 1.985;
double in_xWidth1W = 2.319;
double in_corr1W = 0.123;
double in_weight1 = 0.482 ;
double in_yWidth2N = 1.892;
double in_xWidth2N= 1.886;
double in_corr2N = 0.018;
double in_yWidth2W = 2.141;
double in_xWidth2W = 2.045;
double in_corr2W = 0.212;
double in_weight2 = 0.609;*/

/*double in_yWidth1N= 1.67525;
  double in_xWidth1N= 1.9453;
  double in_corr1N= 0.197133;
  double in_yWidth1W= 2.08182;
  double in_xWidth1W= 2.29246;
  double in_corr1W= 0.299927;
 double in_weight1= 0.826932;
 double in_yWidth2N= 1.99819;
 double in_xWidth2N= 1.8709;
 double in_corr2N= 0.306212;
 double in_yWidth2W= 2.27254;
 double in_xWidth2W= 2.23806;
 double in_corr2W= 0.0156196;
 double in_weight2= 0.264467;*/

     //TFile *f = TFile::Open("./DataAnalysisBunch"+suffix+"DG_with_fineWithCor.root");
 


     TH1F *xwidth1W_h=(TH1F*) f->Get("xwidth1W_h");
     TH1F *xwidth1N_h=(TH1F*) f->Get("xwidth1N_h");
     TH1F *xwidth2W_h=(TH1F*) f->Get("xwidth2W_h");
     TH1F *xwidth2N_h=(TH1F*) f->Get("xwidth2N_h");

     TH1F *ywidth1W_h=(TH1F*) f->Get("ywidth1W_h");
     TH1F *ywidth1N_h=(TH1F*) f->Get("ywidth1N_h");
     TH1F *ywidth2W_h=(TH1F*) f->Get("ywidth2W_h");
     TH1F *ywidth2N_h=(TH1F*) f->Get("ywidth2N_h");

     TH1F *rho1W_h=(TH1F*) f->Get("rho1W_h");
     TH1F *rho1N_h=(TH1F*) f->Get("rho1N_h");
     TH1F *rho2W_h=(TH1F*) f->Get("rho2W_h");
     TH1F *rho2N_h=(TH1F*) f->Get("rho2N_h");

     TH1F *weight1_h=(TH1F*) f->Get("weight1_h");
     TH1F *weight2_h=(TH1F*) f->Get("weight2_h"); 



TH2F *xwidth1W_error_h=(TH2F*) f->Get("xwidth1W_error_h");
     TH2F *xwidth1N_error_h=(TH2F*) f->Get("xwidth1N_error_h");
     TH2F *xwidth2W_error_h=(TH2F*) f->Get("xwidth2W_error_h");
     TH2F *xwidth2N_error_h=(TH2F*) f->Get("xwidth2N_error_h");

     TH2F *ywidth1W_error_h=(TH2F*) f->Get("ywidth1W_error_h");
     TH2F *ywidth1N_error_h=(TH2F*) f->Get("ywidth1N_error_h");
     TH2F *ywidth2W_error_h=(TH2F*) f->Get("ywidth2W_error_h");
     TH2F *ywidth2N_error_h=(TH2F*) f->Get("ywidth2N_error_h");

     TH2F *rho1W_error_h=(TH2F*) f->Get("rho1W_error_h");
     TH2F *rho1N_error_h=(TH2F*) f->Get("rho1N_error_h");
     TH2F *rho2W_error_h=(TH2F*) f->Get("rho2W_error_h");
     TH2F *rho2N_error_h=(TH2F*) f->Get("rho2N_error_h");

     TH2F *weight1_error_h=(TH2F*) f->Get("weight1_error_h");
     TH2F *weight2_error_h=(TH2F*) f->Get("weight2_error_h");
     //TH2F *ywidth1N = (TH2F*) f->Get("y")

     double in_yWidth1N = ywidth1N_h->GetMean();//+r.Uniform(-1.,1.)*ywidth1N_error_h->GetMean();//1.59;//1.8+r.Uniform(0.7);
     double in_xWidth1N = xwidth1N_h->GetMean();//+r.Uniform(-1.,1.)*xwidth1N_error_h->GetMean();//1.78;//1.8+r.Uniform(0.7);
     double in_corr1N = rho1N_h->GetMean();//+r.Uniform(-1.,1.)*rho1N_error_h->GetMean();//-0.115;//-0.4+r.Uniform(0.8);
     double in_yWidth1W = ywidth1W_h->GetMean();//+r.Uniform(-1.,1.)*ywidth1W_error_h->GetMean();//2.39;//1.8+r.Uniform(0.7);
     double in_xWidth1W = xwidth1W_h->GetMean();//+r.Uniform(-1.,1.)*xwidth1W_error_h->GetMean();//3.18;//1.8+r.Uniform(0.7);
     double in_corr1W =rho1W_h->GetMean();//+r.Uniform(-1.,1.)*rho1W_error_h->GetMean(); //-0.23;//-0.4+r.Uniform(0.8);
     double in_weight1 =weight1_h->GetMean();//+r.Uniform(-1.,1.)*weight1_error_h->GetMean(); //0.948;//r.Uniform(1.);

     double in_yWidth2N = ywidth2N_h->GetMean();//+r.Uniform(-1.,1.)*ywidth2N_error_h->GetMean();//1.59;//1.8+r.Uniform(0.7);
     double in_xWidth2N = xwidth2N_h->GetMean();//+r.Uniform(-1.,1.)*xwidth2N_error_h->GetMean();//1.78;//1.8+r.Uniform(0.7);
     double in_corr2N = rho2N_h->GetMean();//+r.Uniform(-1.,1.)*rho2N_error_h->GetMean();//-0.115;//-0.4+r.Uniform(0.8);
     double in_yWidth2W = ywidth2W_h->GetMean();//+r.Uniform(-1.,1.)*ywidth2W_error_h->GetMean();//2.39;//1.8+r.Uniform(0.7);
     double in_xWidth2W = xwidth2W_h->GetMean();//+r.Uniform(-1.,1.)*xwidth2W_error_h->GetMean();//3.18;//1.8+r.Uniform(0.7);
     double in_corr2W =rho2W_h->GetMean();//+r.Uniform(-1.,1.)*rho2W_error_h->GetMean(); //-0.23;//-0.4+r.Uniform(0.8);
     double in_weight2 =weight2_h->GetMean();//+r.Uniform(-1.,1.)*weight2_error_h->GetMean(); //0.948;//r.Uniform(1.);

double in_yWidth1N_error = ywidth1N_error_h->GetMean();//1.59;//1.8+r.Uniform(0.7);
     double in_xWidth1N_error = xwidth1N_error_h->GetMean();//1.78;//1.8+r.Uniform(0.7);
     double in_corr1N_error = rho1N_error_h->GetMean();//-0.115;//-0.4+r.Uniform(0.8);
    double in_yWidth1W_error = ywidth1W_error_h->GetMean();//2.39;//1.8+r.Uniform(0.7);
    double in_xWidth1W_error = xwidth1W_error_h->GetMean();//3.18;//1.8+r.Uniform(0.7);
    double in_corr1W_error =rho1W_error_h->GetMean(); //-0.23;//-0.4+r.Uniform(0.8);
    double in_weight1_error =weight1_error_h->GetMean(); //0.948;//r.Uniform(1.);

double in_yWidth2N_error = ywidth2N_error_h->GetMean();//1.59;//1.8+r.Uniform(0.7);
     double in_xWidth2N_error = xwidth2N_error_h->GetMean();//1.78;//1.8+r.Uniform(0.7);
     double in_corr2N_error = rho2N_error_h->GetMean();//-0.115;//-0.4+r.Uniform(0.8);
    double in_yWidth2W_error = ywidth2W_error_h->GetMean();//2.39;//1.8+r.Uniform(0.7);
    double in_xWidth2W_error = xwidth2W_error_h->GetMean();//3.18;//1.8+r.Uniform(0.7);
    double in_corr2W_error =rho2W_error_h->GetMean(); //-0.23;//-0.4+r.Uniform(0.8);
    double in_weight2_error =weight2_error_h->GetMean(); //0.948;//r.Uniform(1.);

    /*double in_yWidth2N = 1.57;//1.8+r.Uniform(0.7);
    double in_xWidth2N = 1.68;//1.8+r.Uniform(0.7);
    double in_corr2N = -0.313;//-0.4+r.Uniform(0.8);
    double in_yWidth2W =1.54;//1.8+r.Uniform(0.7);
    double in_xWidth2W = 1.7;//1.8+r.Uniform(0.7);
    double in_corr2W = -0.18;//-0.4+r.Uniform(0.8);
    double in_weight2 = 0.288;//r.Uniform(1.);*/
 
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
 
 TH1F* hb2x = new TH1F("hb2x","hb2x",scansteps,-0.5,(scansteps-1)+0.5);
 TH1F* hb2y = new TH1F("hb2y","hb2y",scansteps,-0.5,(scansteps-1)+0.5);
 
 double vtxX = 0;
 double vtxY = 0;
 double vtxXsmeared = 0;
 double vtxYsmeared = 0;
 int scanXn = -100;
 int scanYn = -100;
 
 multBeam->SetNpy(500);
 multBeam->SetNpx(500);

 // Beam 1 fixed, beam 2 scan in x-direction
 for(int scanX = 0; scanX<25; scanX+=1){
   scanXn = scanX;   
   multBeam->SetParameter(9,scanX-12+(scanX*0.000));
   multBeam->SetParameter(1,0.0);
   double integral = multBeam->Integral(-30,30,-30,30)*800;
   double binContent = r.Poisson(integral);
   hb2x->SetBinContent(scanXn+1,binContent);
   cout<<"Beam 1 fixed, beam 2 x-scan, integral is "<<integral<<endl;
 }

 multBeam->SetParameter(9,0.0);
 scanXn = -100;

 // Beam 1 fixed, beam 2 scan in y-direction
 for(int scanY = 0; scanY<25; scanY+=1){
   scanYn = scanY;
   multBeam->SetParameter(10,scanY-12.+(scanY*0.000));
   multBeam->SetParameter(9,0.0);
   double integral = multBeam->Integral(-30,30,-30,30)*800;
   double binContent = r.Poisson(integral);
   hb2y->SetBinContent(scanYn+1,binContent);
   cout<<"Beam 1 fixed, beam 2 y-scan, integral is "<<integral<<endl;
 }

 scanYn = -100;
 multBeam->SetParameter(10,0.0);
 multBeam->SetParameter(9,0.0);
 cout<<"this"<<endl;
  // Fits

  RooRealVar nscan("nscan","nscanpt",0.0,25.);
  nscan.setBins(scansteps);

  // Beam 1 fixed, beam 2 scan in x-direction
  // PDFs and yields
  RooRealVar mu1_b2x("mu1_b2x","mean 1",8.5,14.6);
  //mu1_b2x.setConstant(kTRUE);
  RooRealVar mu2_b2x("mu2_b2x","mean 2",8.2,15.1);
  // mu2_b2x.setConstant(kTRUE);
  RooRealVar sigma1_b2x("sigma1_b2x","sigma 1",1.0,2.9);
  RooRealVar sigma12_b2x("sigma12_b2x","sigma 2",2.1,4.);
  RooGaussian gauss1_b2x("gauss1_b2x","gaussian 1",nscan,mu1_b2x,sigma1_b2x);
  RooGaussian gauss2_b2x("gauss2_b2x","gaussian 2",nscan,mu2_b2x,sigma12_b2x); 
  RooRealVar c("c","c", 0.0, 1.); 
  // RooRealVar gauss1_b2x_area("gauss1_b2x_area","gauss1_b2x_area",0.1,1.);
  // RooRealVar gauss2_b2x_area("gauss2_b2x_area","gauss2_b2x_area",0.1,0.9);
  RooAddPdf dg_b2x("dg_b2x","dg_b2x",RooArgList(gauss1_b2x,gauss2_b2x),c/*RooArgList(gauss1_b2x_area,gauss2_b2x_area)*/);
  // Data
  RooDataHist b2x("b2x","b2x", RooArgSet(nscan),hb2x);
  // Fit
  // RooMsgService::instance().setSilentMode(kTRUE);
  RooFitResult *fitRes_b2x = dg_b2x.fitTo(b2x,Save());

  // Beam 1 fixed, beam 2 scan in x-direction
  // PDFs and yields
  RooRealVar mu1_b2y("mu1_b2y","mean 1",8.1,15.8);
  RooRealVar mu2_b2y("mu2_b2y","mean 2",8.4,15.6);
  RooRealVar sigma1_b2y("sigma1_b2y","sigma 1",1.0,2.5);
  RooRealVar sigma12_b2y("sigma12_b2y","sigma 2",2.11,4.);
  RooGaussian gauss1_b2y("gauss1_b2y","gaussian 1",nscan,mu1_b2y,sigma1_b2y);
  RooGaussian gauss2_b2y("gauss2_b2y","gaussian 2",nscan,mu2_b2y,sigma12_b2y);
  RooRealVar k("k","k", 0.0, 1.); 
  // RooRealVar gauss1_b2y_area("gauss1_b2y_area","gauss1_b2y_area",0.1,1);
  // RooRealVar gauss2_b2y_area("gauss2_b2y_area","gauss2_b2y_area",0.1,0.9);
  RooAddPdf dg_b2y("dg_b2y","dg_b2y",RooArgList(gauss1_b2y,gauss2_b2y),k/*RooArgList(gauss1_b2y_area,gauss2_b2y_area)*/);
  // Data
  RooDataHist b2y("b2y","b2y", RooArgSet(nscan),hb2y);
  // Fit
  RooFitResult *fitRes_b2y = dg_b2y.fitTo(b2y,Save());

  // Calculate overlap diff

  // Beam 1 fixed, beam 2 scan in y-direction
  RooPlot *vtxframe_b2x = nscan.frame();
  vtxframe_b2x->GetYaxis()->SetNdivisions(505);
  vtxframe_b2x->GetYaxis()->SetTitle("# vertices");
  vtxframe_b2x->GetYaxis()->SetTitleOffset(0.5);
  vtxframe_b2x->GetYaxis()->SetTitleSize(0.09);
  vtxframe_b2x->SetTitle("Double Gaussian Fit to Beam 2 X-Scan");
  b2x.plotOn(vtxframe_b2x);
  //dg_b2x.paramOn(vtxframe_b2x);
  dg_b2x.plotOn(vtxframe_b2x,RooFit::Name("b2xfit"));
  double b2x_chiSquared = vtxframe_b2x->chiSquare();
  // Get max val
  RooCurve* fit_b2x = (RooCurve*) vtxframe_b2x->findObject("b2xfit");
  double peak_b2x = findMax(fit_b2x,-0.5,(scansteps-1)+0.5);
  cout << "peak_b2x is " << peak_b2x << endl;
  double area_b2x = fit_b2x->Integral();

  // Beam 1 fixed, beam 2 scan in y-direction
  RooPlot *vtxframe_b2y = nscan.frame();
  vtxframe_b2y->GetYaxis()->SetNdivisions(505);
  vtxframe_b2y->GetYaxis()->SetTitle("# vertices");
  vtxframe_b2y->GetYaxis()->SetTitleOffset(0.5);
  vtxframe_b2y->GetYaxis()->SetTitleSize(0.09);
  vtxframe_b2y->SetTitle("Double Gaussian Fit to Beam 2 Y-Scan");
  b2y.plotOn(vtxframe_b2y);
  //dg_b2y.paramOn(vtxframe_b2y);
  dg_b2y.plotOn(vtxframe_b2y,RooFit::Name("b2yfit"));
  double b2y_chiSquared = vtxframe_b2y->chiSquare();
  // Get max val
  RooCurve* fit_b2y = (RooCurve*) vtxframe_b2y->findObject("b2yfit");
  double peak_b2y = findMax(fit_b2y,-0.5,(scansteps-1)+0.5);
  cout << "peak_b2y is " << peak_b2y << endl;
  double area_b2y = fit_b2y->Integral();

 cout << "True yield = " << hb2x->Integral() << endl;
 double overlapTrue = multBeam->Integral(-30,30,-30,30)/TMath::Power(100,2);
 cout<<"TRUE INTEGRAL "<<overlapTrue<<endl;
 double overlapFit = (peak_b2x/area_b2x)*(peak_b2y/area_b2y);//2*3.1415926*sigma_b2x.getVal()*sigma_b2y.getVal();
 cout<<"FIT INTEGRAL "<<overlapFit<<endl;
 double overlapDiff = (overlapFit - overlapTrue)/overlapTrue;
 cout<<"Overlap DIFF "<<overlapDiff<<endl;
 cout<<"capSigX*capSigY*2pi True "<<(1.*(0.0508*0.0508)/overlapTrue)<<endl;
 cout<<"capSigX*capSigY*2pi Fit "<<(1.*(0.0508*0.0508)/overlapFit)<<endl;

 capSigAeffDiff->Fill((1.*(0.0508*0.0508)/overlapFit));


 cout<<"capSigX :"<<0.0508/(peak_b2x/area_b2x)/2.51<<endl;
 cout<<"capSigY :"<<0.0508/(peak_b2y/area_b2y)/2.51<<endl;
 cout<<"correction factor "<<((1.*(0.00508*0.00508)/overlapTrue)-(1.*(0.00508*0.00508)/overlapFit))/(1.*(0.00508*0.00508)/overlapFit)<<endl;
 cout<<dg_b2y.getCoefNormalization().getRealValue("k")<<endl;

 overDiff->Fill(overlapDiff);

  if(z==1000) {

    TCanvas* cb2x = new TCanvas("cb2x","cb2x");
    drawFit(nscan,vtxframe_b2x,cb2x,"beam2_xscan_TOYS_"+suffix+"_DG.png");
    cb2x->Close();

    TCanvas* cb2y = new TCanvas("cb2y","cb2y");
    drawFit(nscan,vtxframe_b2y,cb2y,"beam2_yscan_TOYS_"+suffix+"_DG.png");
    cb2y->Close();

    // Beam 1 fixed, beam 2 scan in x-direction
    // vtxframe_b2x->Write();
    // Beam 1 fixed, beam 2 scan in y-direction
    // vtxframe_b2y->Write();
  }

 // Save beam parameters and overlap integral calculations
 outTrainTree->Branch("in_yWidth1N",&in_yWidth1N,"in_yWidth1N/D");
 outTrainTree->Branch("in_xWidth1N",&in_xWidth1N,"in_xWidth1N/D");
 outTrainTree->Branch("in_corr1N",&in_corr1N,"in_corr1N/D");
 outTrainTree->Branch("in_yWidth1W",&in_yWidth1W,"in_yWidth1W/D");
 outTrainTree->Branch("in_xWidth1W",&in_xWidth1W,"in_xWidth1W/D");
 outTrainTree->Branch("in_corr1W",&in_corr1W,"in_corr1W/D");
 outTrainTree->Branch("in_weight1",&in_weight1,"in_weight1/D");

 outTrainTree->Branch("in_yWidth2N",&in_yWidth2N,"in_yWidth2N/D");
 outTrainTree->Branch("in_xWidth2N",&in_xWidth2N,"in_xWidth2N/D");
 outTrainTree->Branch("in_corr2N",&in_corr2N,"in_corr2N/D");
 outTrainTree->Branch("in_yWidth2W",&in_yWidth2W,"in_yWidth2W/D");
 outTrainTree->Branch("in_xWidth2W",&in_xWidth2W,"in_xWidth2W/D");
 outTrainTree->Branch("in_corr2W",&in_corr2W,"in_corr2W/D");
 outTrainTree->Branch("in_weight2",&in_weight2,"in_weight2/D");

 double mu1_b2x_fit = mu1_b2x.getVal(), sigma1_b2x_fit = sigma1_b2x.getVal(), mu2_b2x_fit = mu2_b2x.getVal(), sigma2_b2x_fit = sigma12_b2x.getVal();
 double mu1_b2y_fit = mu1_b2y.getVal(), sigma1_b2y_fit = sigma1_b2y.getVal(), mu2_b2y_fit = mu2_b2y.getVal(), sigma2_b2y_fit = sigma12_b2y.getVal();

 outTrainTree->Branch("beam2_xscan_mu1",&mu1_b2x_fit,"beam2_xscan_mu1/D");
 outTrainTree->Branch("beam2_xscan_sigma1",&sigma1_b2x_fit,"beam2_xscan_sigma1/D");
 outTrainTree->Branch("beam2_xscan_mu2",&mu2_b2x_fit,"beam2_xscan_mu2/D");
 outTrainTree->Branch("beam2_xscan_sigma2",&sigma2_b2x_fit,"beam2_xscan_sigma2/D");
 outTrainTree->Branch("beam2_xscan_chiSquared",&b2x_chiSquared,"beam2_xscan_chiSquared/D");

 outTrainTree->Branch("beam2_yscan_mu1",&mu1_b2y_fit,"beam2_yscan_mu1/D");
 outTrainTree->Branch("beam2_yscan_sigma1",&sigma1_b2y_fit,"beam2_yscan_sigma1/D");
 outTrainTree->Branch("beam2_yscan_mu2",&mu2_b2y_fit,"beam2_yscan_mu2/D");
 outTrainTree->Branch("beam2_yscan_sigma2",&sigma2_b2y_fit,"beam2_yscan_sigma2/D");
 outTrainTree->Branch("beam2_yscan_chiSquared",&b2y_chiSquared,"beam2_yscan_chiSquared/D");

 outTrainTree->Branch("overlapTrue",&overlapTrue,"overlapTrue/D");
 outTrainTree->Branch("overlapFit",&overlapFit,"overlapFit/D");
 outTrainTree->Branch("overlapDiff",&overlapDiff,"overlapDiff/D");

 outTrainTree->Fill();
   
 delete hb2x;
 delete hb2y;
   }
 fAna->cd();
 capSigAeffDiff->Write();
 outTrainTree->Write();
 overDiff->Write();
 fAna->Close();	

 //capSigAeffDiff->SaveAs("statErrorVdM.root");											   

}
