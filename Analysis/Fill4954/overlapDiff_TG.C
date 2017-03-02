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

  double beamN1 = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN1,2)))*TMath::Sqrt(TMath::Power(xwidthN1,2.))*TMath::Sqrt(TMath::Power(ywidthN1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN1,2))*(TMath::Power((x[0]-x01)/xwidthN1,2.0)+TMath::Power((x[1]-y01)/ywidthN1,2.0)-2*rhoN1*(x[0]-x01)*(x[1]-y01)/(xwidthN1*ywidthN1)));
  double beamM1 = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoM1,2)))*TMath::Sqrt(TMath::Power(xwidthM1,2.))*TMath::Sqrt(TMath::Power(ywidthM1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoM1,2))*(TMath::Power((x[0]-x01)/xwidthM1,2.0)+TMath::Power((x[1]-y01)/ywidthM1,2.0)-2*rhoM1*(x[0]-x01)*(x[1]-y01)/(xwidthM1*ywidthM1)));
  double beamW1 = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW1,2)))*TMath::Sqrt(TMath::Power(xwidthW1,2.))*TMath::Sqrt(TMath::Power(ywidthW1,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW1,2))*(TMath::Power((x[0]-x01)/xwidthW1,2.0)+TMath::Power((x[1]-y01)/ywidthW1,2.0)-2*rhoW1*(x[0]-x01)*(x[1]-y01)/(xwidthW1*ywidthW1)));
  double beam1 = nw_weight1N * beamN1 + nw_weight1M * beamM1 + nw_weight1W * beamW1;

  double beamN2 = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoN2,2)))*TMath::Sqrt(TMath::Power(xwidthN2,2.))*TMath::Sqrt(TMath::Power(ywidthN2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoN2,2))*(TMath::Power((x[0]-x02)/xwidthN2,2.0)+TMath::Power((x[1]-y02)/ywidthN2,2.0)-2*rhoN2*(x[0]-x02)*(x[1]-y02)/(xwidthN2*ywidthN2)));
  double beamM2 = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoM2,2)))*TMath::Sqrt(TMath::Power(xwidthM2,2.))*TMath::Sqrt(TMath::Power(ywidthM2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoM2,2))*(TMath::Power((x[0]-x02)/xwidthM2,2.0)+TMath::Power((x[1]-y02)/ywidthM2,2.0)-2*rhoM2*(x[0]-x02)*(x[1]-y02)/(xwidthM2*ywidthM2)));
  double beamW2 = 100./(2*pi*TMath::Sqrt((1.-TMath::Power(rhoW2,2)))*TMath::Sqrt(TMath::Power(xwidthW2,2.))*TMath::Sqrt(TMath::Power(ywidthW2,2.)))*TMath::Exp(-0.5/(1.-TMath::Power(rhoW2,2))*(TMath::Power((x[0]-x02)/xwidthW2,2.0)+TMath::Power((x[1]-y02)/ywidthW2,2.0)-2*rhoW2*(x[0]-x02)*(x[1]-y02)/(xwidthW2*ywidthW2)));
  double beam2 = nw_weight2N * beamN2 + nw_weight2M * beamM2 + nw_weight2W * beamW2;

  double product = beam1 * beam2;

  return product;


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

void overlapDiff_TG(TString suffix)
{
    TRandom3 r;
    gStyle->SetOptStat(0);
    r.SetSeed(0);
    int scansteps = 25.;
    double scaling = 0.00458;
    TFile *f = TFile::Open("DataAnalysisBunch"+suffix+"_new_StronRescale.root");

    TH1F *xwidth1W_h=(TH1F*) f->Get("xwidth1W_h");
    TH1F *xwidth1M_h=(TH1F*) f->Get("xwidth1M_h");
    TH1F *xwidth1N_h=(TH1F*) f->Get("xwidth1N_h");
    TH1F *xwidth2W_h=(TH1F*) f->Get("xwidth2W_h");
    TH1F *xwidth2M_h=(TH1F*) f->Get("xwidth2M_h");
    TH1F *xwidth2N_h=(TH1F*) f->Get("xwidth2N_h");

    TH1F *ywidth1W_h=(TH1F*) f->Get("ywidth1W_h");
    TH1F *ywidth1M_h=(TH1F*) f->Get("ywidth1M_h");
    TH1F *ywidth1N_h=(TH1F*) f->Get("ywidth1N_h");
    TH1F *ywidth2W_h=(TH1F*) f->Get("ywidth2W_h");
    TH1F *ywidth2M_h=(TH1F*) f->Get("ywidth2M_h");
    TH1F *ywidth2N_h=(TH1F*) f->Get("ywidth2N_h");

    TH1F *rho1W_h=(TH1F*) f->Get("rho1W_h");
    TH1F *rho1M_h=(TH1F*) f->Get("rho1M_h");
    TH1F *rho1N_h=(TH1F*) f->Get("rho1N_h");
    TH1F *rho2W_h=(TH1F*) f->Get("rho2W_h");
    TH1F *rho2M_h=(TH1F*) f->Get("rho2M_h");
    TH1F *rho2N_h=(TH1F*) f->Get("rho2N_h");

    TH1F *weight1N_h=(TH1F*) f->Get("weight1N_h");
    TH1F *weight2N_h=(TH1F*) f->Get("weight2N_h");
    TH1F *weight1M_h=(TH1F*) f->Get("weight1M_h");
    TH1F *weight2M_h=(TH1F*) f->Get("weight2M_h");

    TH2F *xwidth1W_error_h=(TH2F*) f->Get("xwidth1W_error_h");
    TH2F *xwidth1M_error_h=(TH2F*) f->Get("xwidth1M_error_h");
    TH2F *xwidth1N_error_h=(TH2F*) f->Get("xwidth1N_error_h");
    TH2F *xwidth2W_error_h=(TH2F*) f->Get("xwidth2W_error_h");
    TH2F *xwidth2M_error_h=(TH2F*) f->Get("xwidth2M_error_h");
    TH2F *xwidth2N_error_h=(TH2F*) f->Get("xwidth2N_error_h");

    TH2F *ywidth1W_error_h=(TH2F*) f->Get("ywidth1W_error_h");
    TH2F *ywidth1M_error_h=(TH2F*) f->Get("ywidth1M_error_h");
    TH2F *ywidth1N_error_h=(TH2F*) f->Get("ywidth1N_error_h");
    TH2F *ywidth2W_error_h=(TH2F*) f->Get("ywidth2W_error_h");
    TH2F *ywidth2M_error_h=(TH2F*) f->Get("ywidth2M_error_h");
    TH2F *ywidth2N_error_h=(TH2F*) f->Get("ywidth2N_error_h");

    TH2F *rho1W_error_h=(TH2F*) f->Get("rho1W_error_h");
    TH2F *rho1M_error_h=(TH2F*) f->Get("rho1M_error_h");
    TH2F *rho1N_error_h=(TH2F*) f->Get("rho1N_error_h");
    TH2F *rho2W_error_h=(TH2F*) f->Get("rho2W_error_h");
    TH2F *rho2M_error_h=(TH2F*) f->Get("rho2M_error_h");
    TH2F *rho2N_error_h=(TH2F*) f->Get("rho2N_error_h");

    TH2F *weight1N_error_h=(TH2F*) f->Get("weight1N_error_h");
    TH2F *weight2N_error_h=(TH2F*) f->Get("weight2N_error_h");
    TH2F *weight1M_error_h=(TH2F*) f->Get("weight1M_error_h");
    TH2F *weight2M_error_h=(TH2F*) f->Get("weight2M_error_h");

    double in_yWidth1N = ywidth1N_h->GetMean();
    double in_xWidth1N = xwidth1N_h->GetMean();
    double in_corr1N = rho1N_h->GetMean();
    double in_yWidth1M = ywidth1M_h->GetMean();
    double in_xWidth1M = xwidth1M_h->GetMean();
    double in_corr1M = rho1M_h->GetMean();
    double in_yWidth1W = ywidth1W_h->GetMean();
    double in_xWidth1W = xwidth1W_h->GetMean();
    double in_corr1W =rho1W_h->GetMean();
    double in_weight1N =weight1N_h->GetMean();
    double in_weight1M =weight1M_h->GetMean();

    double in_yWidth2N = ywidth2N_h->GetMean();
    double in_xWidth2N = xwidth2N_h->GetMean();
    double in_corr2N = rho2N_h->GetMean();
    double in_yWidth2M = ywidth2M_h->GetMean();
    double in_xWidth2M = xwidth2M_h->GetMean();
    double in_corr2M = rho2M_h->GetMean();
    double in_yWidth2W = ywidth2W_h->GetMean();
    double in_xWidth2W = xwidth2W_h->GetMean();
    double in_corr2W =rho2W_h->GetMean();
    double in_weight2N =weight2N_h->GetMean();
    double in_weight2M =weight2M_h->GetMean();

    double in_yWidth1N_error = ywidth1N_error_h->GetMean();
    double in_xWidth1N_error = xwidth1N_error_h->GetMean();
    double in_corr1N_error = rho1N_error_h->GetMean();
    double in_yWidth1M_error = ywidth1M_error_h->GetMean();
    double in_xWidth1M_error = xwidth1M_error_h->GetMean();
    double in_corr1M_error = rho1M_error_h->GetMean();
    double in_yWidth1W_error = ywidth1W_error_h->GetMean();
    double in_xWidth1W_error = xwidth1W_error_h->GetMean();
    double in_corr1W_error =rho1W_error_h->GetMean();
    double in_weight1N_error =weight1N_error_h->GetMean();
    double in_weight1M_error =weight1M_error_h->GetMean();

    double in_yWidth2N_error = ywidth2N_error_h->GetMean();
    double in_xWidth2N_error = xwidth2N_error_h->GetMean();
    double in_corr2N_error = rho2N_error_h->GetMean();
    double in_yWidth2M_error = ywidth2M_error_h->GetMean();
    double in_xWidth2M_error = xwidth2M_error_h->GetMean();
    double in_corr2M_error = rho2M_error_h->GetMean();
    double in_yWidth2W_error = ywidth2W_error_h->GetMean();
    double in_xWidth2W_error = xwidth2W_error_h->GetMean();
    double in_corr2W_error =rho2W_error_h->GetMean();
    double in_weight2N_error =weight2N_error_h->GetMean();
    double in_weight2M_error =weight2M_error_h->GetMean();

    TFile *fAna = new TFile("overlapDiff_TOYS_2016_"+suffix+".root","recreate");

    TH1F *overDiff = new TH1F("overDiff","VdM Scan Measured vs True Overlap",200,-0.1,0.1);
    TH1F *capSigAeffDiff = new TH1F("capSigAeffDiff","VdM Scan Measured vs True Overlap",360,0.,0.12);

    TF2 *multBeam = new TF2("multBeam",beamMultTG,-30,30,-30,30,26);

    multBeam->SetParameter(0, 0.0);
    multBeam->SetParameter(1, 0.0);
    multBeam->SetParameter(2,in_xWidth1N);
    multBeam->SetParameter(3,in_yWidth1N);
    multBeam->SetParameter(4,in_corr1N);
    multBeam->SetParameter(5,in_xWidth1M);
    multBeam->SetParameter(6,in_yWidth1M);
    multBeam->SetParameter(7,in_corr1M);
    multBeam->SetParameter(8,in_xWidth1W);
    multBeam->SetParameter(9,in_yWidth1W);
    multBeam->SetParameter(10,in_corr1W);
    multBeam->SetParameter(11,in_weight1N);
    multBeam->SetParameter(12,in_weight1M);

    multBeam->SetParameter(13, 0.0);
    multBeam->SetParameter(14, 0.0);
    multBeam->SetParameter(15,in_xWidth2N);
    multBeam->SetParameter(16,in_yWidth2N);
    multBeam->SetParameter(17,in_corr2N);
    multBeam->SetParameter(18,in_xWidth2M);
    multBeam->SetParameter(19,in_yWidth2M);
    multBeam->SetParameter(20,in_corr2M);
    multBeam->SetParameter(21,in_xWidth2W);
    multBeam->SetParameter(22,in_yWidth2W);
    multBeam->SetParameter(23,in_corr2W);
    multBeam->SetParameter(24,in_weight2N);
    multBeam->SetParameter(25,in_weight2M);

    multBeam->SetNpy(500);
    multBeam->SetNpx(500);

    cout<<"TEST INTEGRAL "<<multBeam->Integral(-30,30,-30,30)<<endl;

    // Save beam parameters and overlap integral calculations
    TTree *outTrainTree = new TTree("outTrainTree","tree of VdM fit parameters");

    outTrainTree->Branch("in_yWidth1N",&in_yWidth1N,"in_yWidth1N/D");
    outTrainTree->Branch("in_xWidth1N",&in_xWidth1N,"in_xWidth1N/D");
    outTrainTree->Branch("in_corr1N",&in_corr1N,"in_corr1N/D");
    outTrainTree->Branch("in_yWidth1M",&in_yWidth1M,"in_yWidth1M/D");
    outTrainTree->Branch("in_xWidth1M",&in_xWidth1M,"in_xWidth1M/D");
    outTrainTree->Branch("in_corr1M",&in_corr1M,"in_corr1M/D");
    outTrainTree->Branch("in_yWidth1W",&in_yWidth1W,"in_yWidth1W/D");
    outTrainTree->Branch("in_xWidth1W",&in_xWidth1W,"in_xWidth1W/D");
    outTrainTree->Branch("in_corr1W",&in_corr1W,"in_corr1W/D");
    outTrainTree->Branch("in_weight1N",&in_weight1N,"in_weight1N/D");
    outTrainTree->Branch("in_weight1M",&in_weight1M,"in_weight1M/D");

    outTrainTree->Branch("in_yWidth2N",&in_yWidth2N,"in_yWidth2N/D");
    outTrainTree->Branch("in_xWidth2N",&in_xWidth2N,"in_xWidth2N/D");
    outTrainTree->Branch("in_corr2N",&in_corr2N,"in_corr2N/D");
    outTrainTree->Branch("in_yWidth2M",&in_yWidth2M,"in_yWidth2M/D");
    outTrainTree->Branch("in_xWidth2M",&in_xWidth2M,"in_xWidth2M/D");
    outTrainTree->Branch("in_corr2M",&in_corr2M,"in_corr2M/D");
    outTrainTree->Branch("in_yWidth2W",&in_yWidth2W,"in_yWidth2W/D");
    outTrainTree->Branch("in_xWidth2W",&in_xWidth2W,"in_xWidth2W/D");
    outTrainTree->Branch("in_corr2W",&in_corr2W,"in_corr2W/D");
    outTrainTree->Branch("in_weight2N",&in_weight2N,"in_weight2N/D");
    outTrainTree->Branch("in_weight2M",&in_weight2M,"in_weight2M/D");

    double mu1_b2x_fit, sigma1_b2x_fit, mu2_b2x_fit, sigma2_b2x_fit;
    double mu1_b2y_fit, sigma1_b2y_fit, mu2_b2y_fit, sigma2_b2y_fit;
    double overlapTrue, overlapFit, overlapDiff, b2x_chiSquared, b2y_chiSquared;

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


    for(int z = 0 ; z < 1000 ; z++){

        TH1F* hb2x = new TH1F("hb2x","hb2x",scansteps,-0.5,(scansteps-1)+0.5);
        TH1F* hb2y = new TH1F("hb2y","hb2y",scansteps,-0.5,(scansteps-1)+0.5);

        // Beam 1 fixed, beam 2 scan in x-direction
        for(int scanX = 0; scanX<25; scanX+=1){
            multBeam->SetParameter(13,scanX-12+(scanX*0.000));
            multBeam->SetParameter(14,0.0);
            double integral = multBeam->Integral(-30,30,-30,30)*800;
            double binContent = r.Poisson(integral);
            hb2x->SetBinContent(scanX+1,binContent);
            cout<<"Beam 1 fixed, beam 2 x-scan, integral is "<<integral<<endl;
        }

        // Beam 1 fixed, beam 2 scan in y-direction
        for(int scanY = 0; scanY<25; scanY+=1){
            multBeam->SetParameter(14,scanY-12.+(scanY*0.000));
            multBeam->SetParameter(13,0.0);
            double integral = multBeam->Integral(-30,30,-30,30)*800;
            double binContent = r.Poisson(integral);
            hb2y->SetBinContent(scanY+1,binContent);
            cout<<"Beam 1 fixed, beam 2 y-scan, integral is "<<integral<<endl;
        }

        multBeam->SetParameter(13,0.0);
        multBeam->SetParameter(14,0.0);

        RooRealVar nscan("nscan","nscanpt",0.0,25.);
        nscan.setBins(scansteps);

        // Beam 1 fixed, beam 2 scan in x-direction
        // PDFs and yields
        RooRealVar mu1_b2x("mu1_b2x","mean 1",8.5,14.6);
        RooRealVar mu2_b2x("mu2_b2x","mean 2",8.2,15.1);
        RooRealVar sigma1_b2x("sigma1_b2x","sigma 1",1.0,2.9);
        RooRealVar sigma12_b2x("sigma12_b2x","sigma 2",2.1,4.);
        RooGaussian gauss1_b2x("gauss1_b2x","gaussian 1",nscan,mu1_b2x,sigma1_b2x);
        RooGaussian gauss2_b2x("gauss2_b2x","gaussian 2",nscan,mu2_b2x,sigma12_b2x);
        RooRealVar c("c","c", 0.0, 1.);
        RooAddPdf dg_b2x("dg_b2x","dg_b2x",RooArgList(gauss1_b2x,gauss2_b2x),c);
        // Data
        RooDataHist b2x("b2x","b2x", RooArgSet(nscan),hb2x);
        // Fit
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
        vtxframe_b2x->SetTitle("Triple Gaussian Fit to Beam 2 X-Scan");
        b2x.plotOn(vtxframe_b2x);
        dg_b2x.plotOn(vtxframe_b2x,RooFit::Name("b2xfit"));
        b2x_chiSquared = vtxframe_b2x->chiSquare();
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
        vtxframe_b2y->SetTitle("Triple Gaussian Fit to Beam 2 Y-Scan");
        b2y.plotOn(vtxframe_b2y);
        dg_b2y.plotOn(vtxframe_b2y,RooFit::Name("b2yfit"));
        b2y_chiSquared = vtxframe_b2y->chiSquare();
        // Get max val
        RooCurve* fit_b2y = (RooCurve*) vtxframe_b2y->findObject("b2yfit");
        double peak_b2y = findMax(fit_b2y,-0.5,(scansteps-1)+0.5);
        cout << "peak_b2y is " << peak_b2y << endl;
        double area_b2y = fit_b2y->Integral();

        cout << "True yield = " << hb2x->Integral() << endl;
        overlapTrue = multBeam->Integral(-30,30,-30,30)/TMath::Power(100,2);
        cout<<"TRUE INTEGRAL "<<overlapTrue<<endl;
        overlapFit = (peak_b2x/area_b2x)*(peak_b2y/area_b2y);
        cout<<"FIT INTEGRAL "<<overlapFit<<endl;
        overlapDiff = (overlapFit - overlapTrue)/overlapTrue;
        cout<<"Overlap DIFF "<<overlapDiff<<endl;
        cout<<"capSigX*capSigY*2pi True "<<(1.*(scaling*scaling)/overlapTrue)<<endl;
        cout<<"capSigX*capSigY*2pi Fit "<<(1.*(scaling*scaling)/overlapFit)<<endl;

        capSigAeffDiff->Fill((1.*(scaling*scaling*scaling*scaling)/overlapFit));

        cout<<"capSigX :"<<scaling/(peak_b2x/area_b2x)/2.51<<endl;
        cout<<"capSigY :"<<scaling/(peak_b2y/area_b2y)/2.51<<endl;
        cout<<"correction factor "<<((1.*(scaling*scaling)/overlapTrue)-(1.*(scaling*scaling)/overlapFit))/(1.*(scaling*scaling)/overlapFit)<<endl;
        cout<<dg_b2y.getCoefNormalization().getRealValue("k")<<endl;

        overDiff->Fill(overlapDiff);

        if(z==999 && false) {
            TCanvas* cb2x = new TCanvas("cb2x","cb2x");
            drawFit(nscan,vtxframe_b2x,cb2x,"beam2_xscan_TOYS_"+suffix+".png");
            cb2x->Close();

            TCanvas* cb2y = new TCanvas("cb2y","cb2y");
            drawFit(nscan,vtxframe_b2y,cb2y,"beam2_yscan_TOYS_"+suffix+".png");
            cb2y->Close();
        }

        mu1_b2x_fit = mu1_b2x.getVal();
        sigma1_b2x_fit = sigma1_b2x.getVal();
        mu2_b2x_fit = mu2_b2x.getVal();
        sigma2_b2x_fit = sigma12_b2x.getVal();
        mu1_b2y_fit = mu1_b2y.getVal();
        sigma1_b2y_fit = sigma1_b2y.getVal();
        mu2_b2y_fit = mu2_b2y.getVal();
        sigma2_b2y_fit = sigma12_b2y.getVal();

        outTrainTree->Fill();

        delete hb2x;
        delete hb2y;
    }
    fAna->cd();
    capSigAeffDiff->Write();
    outTrainTree->Write();
    overDiff->Write();
    fAna->Close();

}
