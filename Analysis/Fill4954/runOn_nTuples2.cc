#include "TH2F.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TMath.h"
#include "TString.h"
//#include "RooFit.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "TRint.h"
//#include "TROOT.h"
//using namespace std;
//#include "RooWorkspace.h"
#include "TTree.h"

using namespace std;

void runOn_nTuples2 (const TString txtFile){

  //41,281,872,1783,2063
  
  TChain chain("lumi/tree");
  
  TString filename;
  ifstream ifs;
  ifs.open(txtFile);
  assert(ifs.is_open());
  string line;
  while (getline(ifs, line)) {
    stringstream ss(line);
    ss >> filename;
    chain.Add(filename);
  }
  ifs.close();
  
  cout <<" Entries: "<<chain.GetEntries() << endl;
  
  int bunchCrossing;
  
  //TFile* file = TFile::Open("root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May182016_274100/ZeroBias2/PCC_VdM_ZeroBias2_274100_ProMay312016_Event_AlwaysTrue/160531_211331/0000/pcc_Data_PixVtx_Event_29.root");
  //TTree *tree = (TTree*) file->Get("lumi/tree");

  UInt_t event;

  int nVtx;
  int vtx_nTrk[200];
  float vtx_x[200]; 
  float vtx_y[200];
  float vtx_z[200];
  float vtx_xError[200]; 
  float vtx_yError[200];
  float vtx_zError[200];
  bool vtx_isGood[200];
  bool vtx_isValid[200];
  bool vtx_isFake[200];
  UInt_t timestamp_begin;
  UInt_t timestamp_end;
 
  chain.SetBranchAddress("bunchCrossing",&bunchCrossing);
  chain.SetBranchAddress("timeStamp_begin",&timestamp_begin);
  chain.SetBranchAddress("timeStamp_end",&timestamp_end);
  chain.SetBranchAddress("event",&event);
  chain.SetBranchAddress("nVtx",&nVtx);
  chain.SetBranchAddress("vtx_x",vtx_x);
  chain.SetBranchAddress("vtx_y",vtx_y);
  chain.SetBranchAddress("vtx_z",vtx_z);
  chain.SetBranchAddress("vtx_xError",vtx_xError);
  chain.SetBranchAddress("vtx_yError",vtx_yError);
  chain.SetBranchAddress("vtx_isGood",vtx_isGood);
  chain.SetBranchAddress("vtx_isValid",vtx_isValid);
  chain.SetBranchAddress("vtx_isFake",vtx_isFake);
  chain.SetBranchAddress("vtx_nTrk",vtx_nTrk);

  TFile *outFile = new TFile("outFile"+txtFile+".root","recreate");

  int begin_X1Movescan = 1464341247;               
  int end_X1Movescan = 1464342191+30;
  int begin_Y1Movescan = 1464342677;
  int end_Y1Movescan = 1464343621+30;
  
  int begin_X2Movescan = 1464344429;
  int end_X2Movescan =  1464345370+30;
  int begin_Y2Movescan = 1464345652;
  int end_Y2Movescan = 1464346604+30;

  int scanX1MoveBegin[19] = {1464341247,1464341300,1464341352, 1464341405, 1464341457, 1464341510, 1464341562, 1464341614, 1464341667, 1464341721, 1464341772, 1464341824, 1464341877, 1464341929, 1464341982, 1464342034, 1464342086, 1464342139, 1464342191};

  int scanX1MoveEnd[19] = {1440472600, 1440472642, 1440472686, 1440472728, 1440472772, 1440472814, 1440472858, 1440472900, 1440472944, 1440472986, 1440473028, 1440473072, 1440473115, 1440473158, 1440473200, 1440473244, 1440473286, 1440473328, 1440473372};
  
  int scanY1MoveBegin[19] = {1464342677, 1464342729, 1464342781, 1464342834, 1464342886, 1464342939, 1464342993, 1464343044, 1464343096, 1464343150, 1464343201,1464343254, 1464343306, 1464343358, 1464343411, 1464343463, 1464343516, 1464343568, 1464343621};

  int scanY1MoveEnd[19] = {1440473682, 1440473727, 1440473770, 1440473815, 1440473857, 1440473901, 1440473945, 1440473988, 1440474032, 1440474076, 1440474119, 1440474163, 1440474207, 1440474250, 1440474294, 1440474339, 1440474382, 1440474427, 1440474470};
  
  int scanX2MoveBegin[19] = {1464344429, 1464344482, 1464344534, 1464344587, 1464344639, 1464344692, 1464344744, 1464344795, 1464344847, 1464344887, 1464344952, 1464345005, 1464345057, 1464345110, 1464345162, 1464345215, 1464345267, 1464345318, 1464345370};

  int scanX2MoveEnd[19] = {1440474803,1440474846,1440474890,1440474934,1440474977,1440475021,1440475063,1440475107,1440475151,1440475195,1440475237,1440475281,1440475324,1440475368,1440475412,1440475455,1440475498,1440475541,1440475585};
  
  int scanY2MoveBegin[19] = {1464345652, 1464345706, 1464345758, 1464345812, 1464345864, 1464345917, 1464345969, 1464346023, 1464346076, 1464346128, 1464346180, 1464346233, 1464346287, 1464346341, 1464346393, 1464346444, 1464346497, 1464346551, 1464346604};

  int scanY2MoveEnd[19] = {1440475894,1440475936,1440475980,1440476022,1440476066,1440476110,1440476152,1440476195,1440476238,1440476281,1440476325,1440476369,1440476411,1440476455,1440476499,1440476541,1440476585,1440476627,1440476670};



  int scanXNr1Begin[25] = {1440463385, 1440463425, 1440463466, 1440463507, 1440463548, 1440463587, 1440463628, 1440463669, 1440463710, 1440463750, 1440463790, 1440463832, 1440463871, 1440463913, 1440463953, 1440463994, 1440464034, 1440464075, 1440464117, 1440464157, 1440464198, 1440464238, 1440464279, 1440464320, 1440464359};
  int scanXNr1End[25] = {1440463412, 1440463453, 1440463494, 1440463535, 1440463576, 1440463616, 1440463657, 1440463698, 1440463737, 1440463778, 1440463819, 1440463860, 1440463900, 1440463941, 1440463982, 1440464023, 1440464063, 1440464104, 1440464145, 1440464186, 1440464227, 1440464268, 1440464307, 1440464348, 1440464388};


int scanYNr1Begin[25] = {1440464725, 1440464764, 1440464805, 1440464846, 1440464887, 1440464927, 1440464968, 1440465009, 1440465051, 1440465091, 1440465133, 1440465172, 1440465213, 1440465254, 1440465295, 1440465335, 1440465376, 1440465417, 1440465458, 1440465499, 1440465539, 1440465580, 1440465622, 1440465663, 1440465704};
int scanYNr1End[25] = {1440464753, 1440464793, 1440464834, 1440464875, 1440464916, 1440464957, 1440464997, 1440465038, 1440465079, 1440465120, 1440465161, 1440465201, 1440465242, 1440465283, 1440465324, 1440465363, 1440465405, 1440465445, 1440465487, 1440465528, 1440465569, 1440465609, 1440465650, 1440465691, 1440465732};



int scanXNr2Begin[25] = {1440467318, 1440467359, 1440467400, 1440467439, 1440467480, 1440467521, 1440467562, 1440467602, 1440467642, 1440467683, 1440467723, 1440467764, 1440467806, 1440467846, 1440467887, 1440467927, 1440467968, 1440468009, 1440468048, 1440468091, 1440468130, 1440468171, 1440468211, 1440468252, 1440468293};
int scanXNr2End[25] = {1440467347, 1440467387, 1440467428, 1440467468, 1440467509, 1440467550, 1440467589, 1440467630, 1440467671, 1440467712, 1440467752, 1440467793, 1440467834, 1440467875, 1440467916, 1440467957, 1440467996, 1440468037, 1440468077, 1440468118, 1440468159, 1440468200, 1440468239, 1440468280, 1440468322};


int scanYNr2Begin[25] = {1440466016, 1440466055, 1440466096, 1440466137, 1440466178, 1440466218, 1440466258, 1440466298, 1440466339, 1440466380, 1440466421, 1440466460, 1440466501, 1440466542, 1440466583, 1440466623, 1440466664, 1440466705, 1440466746, 1440466785, 1440466826, 1440466867, 1440466907, 1440466948, 1440466989};
int scanYNr2End[25] = {1440466043, 1440466084, 1440466125, 1440466166, 1440466205, 1440466246, 1440466287, 1440466328, 1440466368, 1440466408, 1440466449, 1440466489, 1440466530, 1440466571, 1440466612, 1440466652, 1440466693, 1440466733, 1440466773, 1440466814, 1440466855, 1440466896, 1440466935, 1440466976, 1440467017};





 int scanXNr5Begin[25] = {1440478566, 1440478607, 1440478646, 1440478687, 1440478728, 1440478768, 1440478809, 1440478849, 1440478889, 1440478930, 1440478971, 1440479010, 1440479051, 1440479092, 1440479133, 1440479173, 1440479214, 1440479255, 1440479296, 1440479337, 1440479377, 1440479417, 1440479458, 1440479498, 1440479539};
int scanXNr5End[25] = {1440478594, 1440478634, 1440478675, 1440478716, 1440478757, 1440478796, 1440478837, 1440478878, 1440478918, 1440478958, 1440478999, 1440479039, 1440479080, 1440479121, 1440479162, 1440479203, 1440479243, 1440479283, 1440479324, 1440479364, 1440479405, 1440479446, 1440479487, 1440479527, 1440479567};


int scanYNr5Begin[25] = {1440479880, 1440479921, 1440479962, 1440480002, 1440480042, 1440480083, 1440480123, 1440480164, 1440480205, 1440480244, 1440480285, 1440480326, 1440480367, 1440480406, 1440480447, 1440480488, 1440480530, 1440480571, 1440480612, 1440480652, 1440480693, 1440480734, 1440480775, 1440480817, 1440480856};
int scanYNr5End[25] = {1440479909, 1440479950, 1440479989, 1440480030, 1440480071, 1440480112, 1440480151, 1440480192, 1440480233, 1440480273, 1440480314, 1440480355, 1440480396, 1440480435, 1440480476, 1440480517, 1440480558, 1440480598, 1440480639, 1440480681, 1440480722, 1440480763, 1440480804, 1440480845, 1440480885};




  int LSC_X1_begin[11]={1440475853,1440476041,1440476165,1440476280,1440476402,1440476519,1440476654,1440476784,1440476909,1440477027,1440477158};

  int LSC_X1_end[11]={1440475895,1440476102,1440476226,1440476344,1440476460,1440476603,1440476730,1440476855,1440476976,1440477116,1440477316};

  int LSC_Y1_begin[11]={1440477575,1440477719,1440477831,1440477952,1440478069,1440478190,1440478312,1440478429,1440478550,1440478671,1440478790};

  int LSC_Y1_end[11]={1440477652,1440477786,1440477904,1440478024,1440478145,1440478263,1440478384,1440478505,1440478624,1440478745,1440478863};




  TString scanN[25] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"};

  TH1F *xError_hist = new TH1F("xError_hist","xError_hist",1000,0.,0.1);
  TH1F *yError_hist = new TH1F("yError_hist","yError_hist",1000,0.,0.1);
  

  int nbinsxy = 10.*19.;
  int nbinsxy2 = 10.*19.;

  TH2F *Beam2MoveX_bunch41Add = new TH2F("Beam2MoveX_bunch41Add","Beam2MoveX_bunch41Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam2MoveY_bunch41Add = new TH2F("Beam2MoveY_bunch41Add","Beam2MoveY_bunch41Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveX_bunch41Add = new TH2F("Beam1MoveX_bunch41Add","Beam1MoveX_bunch41Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveY_bunch41Add = new TH2F("Beam1MoveY_bunch41Add","Beam1MoveY_bunch41Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);

  TH2F *Beam2MoveX_bunch281Add = new TH2F("Beam2MoveX_bunch281Add","Beam2MoveX_bunch281Add",nbinsxy,-10.,10., nbinsxy,-10., 10.);
  TH2F *Beam2MoveY_bunch281Add = new TH2F("Beam2MoveY_bunch281Add","Beam2MoveY_bunch281Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveX_bunch281Add = new TH2F("Beam1MoveX_bunch281Add","Beam1MoveX_bunch281Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveY_bunch281Add = new TH2F("Beam1MoveY_bunch281Add","Beam1MoveY_bunch281Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);

  TH2F *Beam2MoveX_bunch872Add = new TH2F("Beam2MoveX_bunch872Add","Beam2MoveX_bunch872Add",nbinsxy,-10.,10., nbinsxy,-10., 10.);
  TH2F *Beam2MoveY_bunch872Add = new TH2F("Beam2MoveY_bunch872Add","Beam2MoveY_bunch872Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveX_bunch872Add = new TH2F("Beam1MoveX_bunch872Add","Beam1MoveX_bunch872Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveY_bunch872Add = new TH2F("Beam1MoveY_bunch872Add","Beam1MoveY_bunch872Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);


  TH2F *Beam2MoveX_bunch1783Add = new TH2F("Beam2MoveX_bunch1783Add","Beam2MoveX_bunch1783Add",nbinsxy,-10.,10., nbinsxy,-10., 10.);
  TH2F *Beam2MoveY_bunch1783Add = new TH2F("Beam2MoveY_bunch1783Add","Beam2MoveY_bunch1783Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveX_bunch1783Add = new TH2F("Beam1MoveX_bunch1783Add","Beam1MoveX_bunch1783Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveY_bunch1783Add = new TH2F("Beam1MoveY_bunch1783Add","Beam1MoveY_bunch1783Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);

  TH2F *Beam2MoveX_bunch2063Add = new TH2F("Beam2MoveX_bunch2063Add","Beam2MoveX_bunch2063Add",nbinsxy,-10.,10., nbinsxy,-10., 10.);
  TH2F *Beam2MoveY_bunch2063Add = new TH2F("Beam2MoveY_bunch2063Add","Beam2MoveY_bunch2063Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveX_bunch2063Add = new TH2F("Beam1MoveX_bunch2063Add","Beam1MoveX_bunch2063Add",nbinsxy,-10.,10., nbinsxy,-10.,10.);
  TH2F *Beam1MoveY_bunch2063Add = new TH2F("Beam1MoveY_bunch2063Add","Beam1MoveY_bunch2063Add",nbinsxy,-10.,10., nbinsxy,-10.,10. );

  TH2F *Vtx_41_zy = new TH2F("Vtx_41_zy","Vtx_41zy",nbinsxy*4.,0.,0.2, nbinsxy*4.,-20., 20.);	
  TH2F *Vtx_281_zy = new TH2F("Vtx_281_zy","Vtx_281zy",nbinsxy*4.,0.,0.2, nbinsxy*4.,-20., 20.);	
  TH2F *Vtx_872_zy = new TH2F("Vtx_872_zy","Vtx_872zy",nbinsxy*4.,0.,0.2, nbinsxy*4.,-20., 20.);
  TH2F *Vtx_1783_zy = new TH2F("Vtx_1783_zy","Vtx_1783zy",nbinsxy*4.,0.,0.2, nbinsxy*4.,-20., 20.);	
  TH2F *Vtx_2063_zy = new TH2F("Vtx_2063_zy","Vtx_2063zy",nbinsxy*4.,0.,0.2, nbinsxy*4.,-20., 20.);

  TH2F *scanX1Move[19];
  TH2F *scanY1Move[19];
  TH2F *scanX2Move[19];
  TH2F *scanY2Move[19];

  TH2F *scanX1Move_b41[19];
  TH2F *scanY1Move_b41[19];
  TH2F *scanX2Move_b41[19];
  TH2F *scanY2Move_b41[19];

  TH2F *scanXNr1_b41[25];
  TH2F *scanYNr1_b41[25];
  TH2F *scanXNr2_b41[25];
  TH2F *scanYNr2_b41[25];
  TH2F *scanXNr5_b41[25];
  TH2F *scanYNr5_b41[25];
  cout <<" Entries: "<<chain.GetEntries() << endl;
  TH2F *scanX1Move_b281[19];
  TH2F *scanY1Move_b281[19];
  TH2F *scanX2Move_b281[19];
  TH2F *scanY2Move_b281[19];

  TH2F *scanXNr1_b281[25];
  TH2F *scanYNr1_b281[25];
  TH2F *scanXNr2_b281[25];
  TH2F *scanYNr2_b281[25];
  TH2F *scanXNr5_b281[25];
  TH2F *scanYNr5_b281[25];
  cout <<" Entries: "<<chain.GetEntries() << endl;

  TH2F *scanX1Move_b872[19];
  TH2F *scanY1Move_b872[19];
  TH2F *scanX2Move_b872[19];
  TH2F *scanY2Move_b872[19];

  TH2F *scanXNr1_b872[25];
  TH2F *scanYNr1_b872[25];
  TH2F *scanXNr2_b872[25];
  TH2F *scanYNr2_b872[25];
  TH2F *scanXNr5_b872[25];
  TH2F *scanYNr5_b872[25];


  TH2F *scanX1Move_b1783[19];
  TH2F *scanY1Move_b1783[19];
  TH2F *scanX2Move_b1783[19];
  TH2F *scanY2Move_b1783[19];

  TH2F *scanXNr1_b1783[25];
  TH2F *scanYNr1_b1783[25];
  TH2F *scanXNr2_b1783[25];
  TH2F *scanYNr2_b1783[25];
  TH2F *scanXNr5_b1783[25];
  TH2F *scanYNr5_b1783[25];

  TH2F *scanX1Move_b2063[19];
  TH2F *scanY1Move_b2063[19];
  TH2F *scanX2Move_b2063[19];
  TH2F *scanY2Move_b2063[19];

  TH2F *scanXNr1_b2063[25];
  TH2F *scanYNr1_b2063[25];
  TH2F *scanXNr2_b2063[25];
  TH2F *scanYNr2_b2063[25];
  TH2F *scanXNr5_b2063[25];
  TH2F *scanYNr5_b2063[25];



  TH2F *lscX1[11];
  TH2F *lscY1[11];
  cout <<" Entries: "<<chain.GetEntries() << endl;
  for(int i = 0; i<25;i++){

   
    scanXNr1_b41[i] = new TH2F("scanXNr1_b41_"+scanN[i],"scanXNr1_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
  cout <<" Entries: hall"<<chain.GetEntries() << endl;
    scanXNr1_b41[i]->SetOption("colz");
  cout <<" Entries: hall"<<chain.GetEntries() << endl;
    scanXNr1_b281[i] = new TH2F("scanXNr1_b281_"+scanN[i],"scanXNr1_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr1_b281[i]->SetOption("colz");
    scanXNr1_b872[i] = new TH2F("scanXNr1_b872_"+scanN[i],"scanXNr1_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr1_b872[i]->SetOption("colz");
    scanXNr1_b1783[i] = new TH2F("scanXNr1_b1783_"+scanN[i],"scanXNr1_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr1_b1783[i]->SetOption("colz");
    scanXNr1_b2063[i] = new TH2F("scanXNr1_b2063_"+scanN[i],"scanXNr1_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);     
scanXNr1_b2063[i]->SetOption("colz");
  cout <<" Entries: hall"<<chain.GetEntries() << endl;
    scanYNr1_b41[i] = new TH2F("scanYNr1_b41_"+scanN[i],"scanYNr1_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr1_b41[i]->SetOption("colz");
    scanYNr1_b281[i] = new TH2F("scanYNr1_b281_"+scanN[i],"scanYNr1_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr1_b281[i]->SetOption("colz");
    scanYNr1_b872[i] = new TH2F("scanYNr1_b872_"+scanN[i],"scanYNr1_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr1_b872[i]->SetOption("colz");
    scanYNr1_b1783[i] = new TH2F("scanYNr1_b1783_"+scanN[i],"scanYNr1_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr1_b1783[i]->SetOption("colz");
    scanYNr1_b2063[i] = new TH2F("scanYNr1_b2063_"+scanN[i],"scanYNr1_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);     
scanYNr1_b2063[i]->SetOption("colz");
  cout <<" Entries: hall"<<chain.GetEntries() << endl;
    scanXNr2_b41[i] = new TH2F("scanXNr2_b41_"+scanN[i],"scanXNr2_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr2_b41[i]->SetOption("colz");
    scanXNr2_b281[i] = new TH2F("scanXNr2_b281_"+scanN[i],"scanXNr2_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr2_b281[i]->SetOption("colz");
    scanXNr2_b872[i] = new TH2F("scanXNr2_b872_"+scanN[i],"scanXNr2_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr2_b872[i]->SetOption("colz");
    scanXNr2_b1783[i] = new TH2F("scanXNr2_b1783_"+scanN[i],"scanXNr2_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr2_b1783[i]->SetOption("colz");
    scanXNr2_b2063[i] = new TH2F("scanXNr2_b2063_"+scanN[i],"scanXNr2_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);     
scanXNr2_b2063[i]->SetOption("colz");

    scanYNr2_b41[i] = new TH2F("scanYNr2_b41_"+scanN[i],"scanYNr2_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr2_b41[i]->SetOption("colz");
    scanYNr2_b281[i] = new TH2F("scanYNr2_b281_"+scanN[i],"scanYNr2_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr2_b281[i]->SetOption("colz");
    scanYNr2_b872[i] = new TH2F("scanYNr2_b872_"+scanN[i],"scanYNr2_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr2_b872[i]->SetOption("colz");
    scanYNr2_b1783[i] = new TH2F("scanYNr2_b1783_"+scanN[i],"scanYNr2_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr2_b1783[i]->SetOption("colz");
    scanYNr2_b2063[i] = new TH2F("scanYNr2_b2063_"+scanN[i],"scanYNr2_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.); 
    scanYNr2_b2063[i]->SetOption("colz");

    scanXNr5_b41[i] = new TH2F("scanXNr5_b41_"+scanN[i],"scanXNr5_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr5_b41[i]->SetOption("colz");
    scanXNr5_b281[i] = new TH2F("scanXNr5_b281_"+scanN[i],"scanXNr5_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr5_b281[i]->SetOption("colz");
    scanXNr5_b872[i] = new TH2F("scanXNr5_b872_"+scanN[i],"scanXNr5_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr5_b872[i]->SetOption("colz");
    scanXNr5_b1783[i] = new TH2F("scanXNr5_b1783_"+scanN[i],"scanXNr5_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanXNr5_b1783[i]->SetOption("colz");
    scanXNr5_b2063[i] = new TH2F("scanXNr5_b2063_"+scanN[i],"scanXNr5_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);    

    scanYNr5_b41[i] = new TH2F("scanYNr5_b41_"+scanN[i],"scanYNr5_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr5_b41[i]->SetOption("colz");
    scanYNr5_b281[i] = new TH2F("scanYNr5_b281_"+scanN[i],"scanYNr5_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr5_b281[i]->SetOption("colz");
    scanYNr5_b872[i] = new TH2F("scanYNr5_b872_"+scanN[i],"scanYNr5_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr5_b872[i]->SetOption("colz");
    scanYNr5_b1783[i] = new TH2F("scanYNr5_b1783_"+scanN[i],"scanYNr5_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanYNr5_b1783[i]->SetOption("colz");
    scanYNr5_b2063[i] = new TH2F("scanYNr5_b2063_"+scanN[i],"scanYNr5_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);    
 scanYNr5_b2063[i]->SetOption("colz");
  //cout <<" Entries: "<<chain.GetEntries() << endl;
    if(i<19){

    scanX1Move[i] = new TH2F("scanX1Move_"+scanN[i],"scanX1Move_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX1Move[i]->SetOption("colz");
    scanX1Move_b41[i] = new TH2F("scanX1Move_b41_"+scanN[i],"scanX1Move_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX1Move_b41[i]->SetOption("colz");
    scanX1Move_b281[i] = new TH2F("scanX1Move_b281_"+scanN[i],"scanX1Move_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX1Move_b281[i]->SetOption("colz");
    scanX1Move_b872[i] = new TH2F("scanX1Move_b872_"+scanN[i],"scanX1Move_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX1Move_b872[i]->SetOption("colz");
    scanX1Move_b1783[i] = new TH2F("scanX1Move_b1783_"+scanN[i],"scanX1Move_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX1Move_b1783[i]->SetOption("colz");
    scanX1Move_b2063[i] = new TH2F("scanX1Move_b2063_"+scanN[i],"scanX1Move_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX1Move_b2063[i]->SetOption("colz");
    
    scanY1Move[i] = new TH2F("scanY1Move_"+scanN[i],"scanY1Move_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY1Move[i]->SetOption("colz");
    scanY1Move_b41[i] = new TH2F("scanY1Move_b41_"+scanN[i],"scanY1Move_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY1Move_b41[i]->SetOption("colz");
    scanY1Move_b281[i] = new TH2F("scanY1Move_b281_"+scanN[i],"scanY1Move_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY1Move_b281[i]->SetOption("colz");
    scanY1Move_b872[i] = new TH2F("scanY1Move_b872_"+scanN[i],"scanY1Move_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY1Move_b872[i]->SetOption("colz");
    scanY1Move_b1783[i] = new TH2F("scanY1Move_b1783_"+scanN[i],"scanY1Move_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY1Move_b1783[i]->SetOption("colz");
    scanY1Move_b2063[i] = new TH2F("scanY1Move_b2063_"+scanN[i],"scanY1Move_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY1Move_b2063[i]->SetOption("colz");

    scanX2Move[i] = new TH2F("scanX2Move_"+scanN[i],"scanX2Move_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX2Move[i]->SetOption("colz");
    scanX2Move_b41[i] = new TH2F("scanX2Move_b41_"+scanN[i],"scanX2Move_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX2Move_b41[i]->SetOption("colz");
    scanX2Move_b281[i] = new TH2F("scanX2Move_b281_"+scanN[i],"scanX2Move_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX2Move_b281[i]->SetOption("colz");
    scanX2Move_b872[i] = new TH2F("scanX2Move_b872_"+scanN[i],"scanX2Move_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX2Move_b872[i]->SetOption("colz");
    scanX2Move_b1783[i] = new TH2F("scanX2Move_b1783_"+scanN[i],"scanX2Move_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX2Move_b1783[i]->SetOption("colz");
    scanX2Move_b2063[i] = new TH2F("scanX2Move_b2063_"+scanN[i],"scanX2Move_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanX2Move_b2063[i]->SetOption("colz");

    scanY2Move[i] = new TH2F("scanY2Move_"+scanN[i],"scanY2Move_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY2Move[i]->SetOption("colz");
    scanY2Move_b41[i] = new TH2F("scanY2Move_b41_"+scanN[i],"scanY2Move_b41_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY2Move_b41[i]->SetOption("colz");
    scanY2Move_b281[i] = new TH2F("scanY2Move_b281_"+scanN[i],"scanY2Move_b281_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY2Move_b281[i]->SetOption("colz");
    scanY2Move_b872[i] = new TH2F("scanY2Move_b872_"+scanN[i],"scanY2Move_b872_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY2Move_b872[i]->SetOption("colz");
    scanY2Move_b1783[i] = new TH2F("scanY2Move_b1783_"+scanN[i],"scanY2Move_b1783_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY2Move_b1783[i]->SetOption("colz");
    scanY2Move_b2063[i] = new TH2F("scanY2Move_b2063_"+scanN[i],"scanY2Move_b2063_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    scanY2Move_b2063[i]->SetOption("colz");

    }
    
  }

for(int i = 0; i<11;i++){
    lscX1[i] = new TH2F("lscX1_"+scanN[i],"lscX1_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    lscX1[i]->SetOption("colz");
    
    lscY1[i] = new TH2F("lscY1_"+scanN[i],"lscY1_"+scanN[i],nbinsxy2,-10.,10.,nbinsxy2,-10.,10.);
    lscY1[i]->SetOption("colz");
  }

  

  TH2F *scanX1MoveProfile = new TH2F("scanX1MoveProfile","scanX1MoveProfile",19,0.,19.,200,-10.*0.0046+0.1,10.*0.0046+0.1);

  TH2F *scanY1MoveProfile = new TH2F("scanY1MoveProfile","scanY1MoveProfile",19,0.,19.,200,-10.*0.0046+0.1,10.*0.0046+0.1);

  TH2F *scanX2MoveProfile = new TH2F("scanX2MoveProfile","scanX2MoveProfile",19,0.,19.,200,-10.*0.0046+0.1,10.*0.0046+0.1);

  TH2F *scanY2MoveProfile = new TH2F("scanY2MoveProfile","scanY2MoveProfile",19,0.,19.,200,-10.*0.0046+0.1,10.*0.0046+0.1);


  for(int ievent = 0; ievent < chain.GetEntries(); ievent++){
    
    if(ievent % 10000 == 0)
      cout<<ievent<<endl;
 //cout<<"in"<<endl;
    chain.GetEntry(ievent);
   

    //cout<<bunchCrossing<<endl;

    ////cout<<timestamp_begin<<endl;
    ////cout<<timestamp_end<<endl<<endl;
    
   
    if(nVtx>0)
      {
      
      for(int ivtx = 0; ivtx<nVtx; ivtx++)
      	{

	 //std:://cout<<"Hi_1"<<std::endl;
	  if(vtx_nTrk[ivtx]>14 && vtx_isGood[ivtx] && !vtx_isFake[ivtx] /*&& (vtx_xError[ivtx] > 0.0047 && vtx_yError[ivtx] > 0.0047)*/){
	     //cout<<"Hi_2222222222"<<std::endl;
	     //cout<<timestamp_begin<<std::endl;
	    xError_hist->Fill(vtx_xError[ivtx]);
	    yError_hist->Fill(vtx_yError[ivtx]);

	    float xVtx = (vtx_x[ivtx] - 0.073 + 0.0077)/0.00608 ;//+ 1.52;
	    float yVtx = (vtx_y[ivtx] - 0.0985 + 0.0031)/0.00608 ;//+ 0.53;



	if(bunchCrossing == 41)
	  {Vtx_41_zy->Fill(vtx_y[ivtx],vtx_z[ivtx]);}
	if(bunchCrossing == 281)
	  {Vtx_281_zy->Fill(vtx_y[ivtx],vtx_z[ivtx]);}
	if(bunchCrossing == 872)
	  {Vtx_872_zy->Fill(vtx_y[ivtx],vtx_z[ivtx]);}
	if(bunchCrossing == 1783)
	  {Vtx_1783_zy->Fill(vtx_y[ivtx],vtx_z[ivtx]);}
	if(bunchCrossing == 2063)
	  {Vtx_2063_zy->Fill(vtx_y[ivtx],vtx_z[ivtx]);}

	  

	    for(int n = 0; n<25 ; n++)
	      {

		if(n<19.){ 

		  if((timestamp_begin >= scanX1MoveBegin[n]) && (timestamp_begin<=scanX1MoveBegin[n]+25 /*scanX1MoveEnd[n]*/))
		  {
		    	    //cout<<"hallo1"<<std::endl;
		    scanX1Move[n]->Fill(xVtx,yVtx);
		    scanX1MoveProfile->Fill((float)n,vtx_y[ivtx]);


		   if(bunchCrossing == 41)
	             {scanX1Move_b41[n]->Fill(xVtx,yVtx);
		      Beam1MoveX_bunch41Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanX1Move_b281[n]->Fill(xVtx,yVtx);
		      Beam1MoveX_bunch281Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanX1Move_b872[n]->Fill(xVtx,yVtx);
		       Beam1MoveX_bunch872Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanX1Move_b1783[n]->Fill(xVtx,yVtx);
		       Beam1MoveX_bunch1783Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanX1Move_b2063[n]->Fill(xVtx,yVtx);
		       Beam1MoveX_bunch2063Add->Fill(xVtx,yVtx);}
		  }

		  if((timestamp_begin >= scanY1MoveBegin[n]) && (timestamp_begin<= scanY1MoveBegin[n]+25/*scanY1MoveEnd[n]*/))
		  {

		    	 //cout<<"hallo2"<<std::endl;
		    scanY1Move[n]->Fill(xVtx,yVtx);
		    scanY1MoveProfile->Fill((float)n,vtx_x[ivtx]);


		   if(bunchCrossing == 41)
	             {scanY1Move_b41[n]->Fill(xVtx,yVtx);
		      Beam1MoveY_bunch41Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanY1Move_b281[n]->Fill(xVtx,yVtx);
		      Beam1MoveY_bunch281Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanY1Move_b872[n]->Fill(xVtx,yVtx);
		      Beam1MoveY_bunch872Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanY1Move_b1783[n]->Fill(xVtx,yVtx);
		      Beam1MoveY_bunch1783Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanY1Move_b2063[n]->Fill(xVtx,yVtx);
		      Beam1MoveY_bunch2063Add->Fill(xVtx,yVtx);}
		  }

		if((timestamp_begin >= scanX2MoveBegin[n]) && (timestamp_begin<= scanX2MoveBegin[n]+25))
		  {

		    	 //cout<<"hallo3"<<std::endl;
		    scanX2Move[n]->Fill(xVtx,yVtx);
		    scanX2MoveProfile->Fill((float)n,vtx_y[ivtx]);


		   if(bunchCrossing == 41)
	             {scanX2Move_b41[n]->Fill(xVtx,yVtx);
		      Beam2MoveX_bunch41Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanX2Move_b281[n]->Fill(xVtx,yVtx);
		       Beam2MoveX_bunch281Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanX2Move_b872[n]->Fill(xVtx,yVtx);
		       Beam2MoveX_bunch872Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanX2Move_b1783[n]->Fill(xVtx,yVtx);
		       Beam2MoveX_bunch1783Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanX2Move_b2063[n]->Fill(xVtx,yVtx);
		       Beam2MoveX_bunch2063Add->Fill(xVtx,yVtx);}
		  }

		if((timestamp_begin >= scanY2MoveBegin[n]) && (timestamp_begin<= scanY2MoveBegin[n]+25))
		  {

		    	 //cout<<"hallo4"<<std::endl;
		    scanY2Move[n]->Fill(xVtx,yVtx);
		    scanY2MoveProfile->Fill((float)n,vtx_x[ivtx]);

	           if(bunchCrossing == 41)
	             {scanY2Move_b41[n]->Fill(xVtx,yVtx);
		      Beam2MoveY_bunch41Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanY2Move_b281[n]->Fill(xVtx,yVtx);
		      Beam2MoveY_bunch281Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanY2Move_b872[n]->Fill(xVtx,yVtx);
		      Beam2MoveY_bunch872Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanY2Move_b1783[n]->Fill(xVtx,yVtx);
		      Beam2MoveY_bunch1783Add->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanY2Move_b2063[n]->Fill(xVtx,yVtx);
		      Beam2MoveY_bunch2063Add->Fill(xVtx,yVtx);}

		    
		  }

		}

		 if((timestamp_begin >= scanXNr1Begin[n]) && (timestamp_begin<= scanXNr1End[n]))
		  {

	           if(bunchCrossing == 41)
	             {scanXNr1_b41[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanXNr1_b281[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanXNr1_b872[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanXNr1_b1783[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanXNr1_b2063[n]->Fill(xVtx,yVtx);}
		  }
            
	         if((timestamp_begin >= scanYNr1Begin[n]) && (timestamp_begin<= scanYNr1End[n]))
		  {

	           if(bunchCrossing == 41)
	             {scanYNr1_b41[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanYNr1_b281[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanYNr1_b872[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanYNr1_b1783[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanYNr1_b2063[n]->Fill(xVtx,yVtx);}
		  }
            

		if((timestamp_begin >= scanXNr2Begin[n]) && (timestamp_begin<= scanXNr2End[n]))
		  {

	           if(bunchCrossing == 41)
	             {scanXNr2_b41[n]->Fill(xVtx,yVtx);
       		      }
	           if(bunchCrossing == 281)
	             {scanXNr2_b281[n]->Fill(xVtx,yVtx);
		      }
	           if(bunchCrossing == 872)
	             {scanXNr2_b872[n]->Fill(xVtx,yVtx);
		      }
	           if(bunchCrossing == 1783)
	             {scanXNr2_b1783[n]->Fill(xVtx,yVtx);
		      }
	           if(bunchCrossing == 2063)
	             {scanXNr2_b2063[n]->Fill(xVtx,yVtx);
		       }
		  }
            
	         if((timestamp_begin >= scanYNr2Begin[n]) && (timestamp_begin<= scanYNr2End[n]))
		  {

	           if(bunchCrossing == 41)
	             {scanYNr2_b41[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanYNr2_b281[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanYNr2_b872[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanYNr2_b1783[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanYNr2_b2063[n]->Fill(xVtx,yVtx);}
		  }


		if((timestamp_begin >= scanXNr5Begin[n]) && (timestamp_begin<= scanXNr5End[n]))
		  {

	           if(bunchCrossing == 41)
	             {scanXNr5_b41[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanXNr5_b281[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanXNr5_b872[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanXNr5_b1783[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanXNr5_b2063[n]->Fill(xVtx,yVtx);}
		  }
            
	         if((timestamp_begin >= scanYNr5Begin[n]) && (timestamp_begin<= scanYNr5End[n]))
		  {

	           if(bunchCrossing == 41)
	             {scanYNr5_b41[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 281)
	             {scanYNr5_b281[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 872)
	             {scanYNr5_b872[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 1783)
	             {scanYNr5_b1783[n]->Fill(xVtx,yVtx);}
	           if(bunchCrossing == 2063)
	             {scanYNr5_b2063[n]->Fill(xVtx,yVtx);}
		  }

	    

	      }








	    for(int k=0;k<11;k++)
	      {
		
	      if(timestamp_begin>=LSC_X1_begin[k] && timestamp_begin<=LSC_X1_end[k])
		{
		  lscX1[k]->Fill(xVtx,yVtx);
		}

	      if(timestamp_begin>=LSC_Y1_begin[k] && timestamp_begin<=LSC_Y1_end[k])
		{
		  lscY1[k]->Fill(xVtx,yVtx);
		}

	      }

	    

	    /*


	    if(timestamp_begin> begin_X1Movescan &&  timestamp_begin<end_X1Movescan)
	      {
		if(bunchCrossing == 41)
		  Beam1MoveX_bunch41Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 281)
		  Beam1MoveX_bunch281Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 872)
		  Beam1MoveX_bunch872Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 1783)
		  Beam1MoveX_bunch1783Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 2063)
		  Beam1MoveX_bunch2063Add->Fill(xVtx,yVtx);		
	      }


	    if(timestamp_begin> begin_Y1Movescan &&  timestamp_begin<end_Y1Movescan)
	      {
		if(bunchCrossing == 41)
		  Beam1MoveY_bunch41Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 281)
		  Beam1MoveY_bunch281Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 872)
		  Beam1MoveY_bunch872Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 1783)
		  Beam1MoveY_bunch1783Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 2063)
		  Beam1MoveY_bunch2063Add->Fill(xVtx,yVtx);		
	      }

	    if(timestamp_begin> begin_X2Movescan &&  timestamp_begin<end_X2Movescan)
	      {
		if(bunchCrossing == 41)
		  Beam2MoveX_bunch41Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 281)
		  Beam2MoveX_bunch281Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 872)
		  Beam2MoveX_bunch872Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 1783)
		  Beam2MoveX_bunch1783Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 2063)
		  Beam2MoveX_bunch2063Add->Fill(xVtx,yVtx);		
	      }


	    if(timestamp_begin> begin_Y2Movescan &&  timestamp_begin<end_Y2Movescan)
	      {
		if(bunchCrossing == 41)
		  Beam2MoveY_bunch41Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 281)
		  Beam2MoveY_bunch281Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 872)
		  Beam2MoveY_bunch872Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 1783)
		  Beam2MoveY_bunch1783Add->Fill(xVtx,yVtx);
		if(bunchCrossing == 2063)
		  Beam2MoveY_bunch2063Add->Fill(xVtx,yVtx);		
	      }

	    */
	   
	    
	
	
	  }
	}
    
      

    }
  }


  Beam2MoveX_bunch41Add->SetOption("colz");
  Beam2MoveY_bunch41Add->SetOption("colz");
  Beam1MoveX_bunch41Add->SetOption("colz");
  Beam1MoveY_bunch41Add->SetOption("colz");

  Beam2MoveX_bunch281Add->SetOption("colz");
  Beam2MoveY_bunch281Add->SetOption("colz");
  Beam1MoveX_bunch281Add->SetOption("colz");
  Beam1MoveY_bunch281Add->SetOption("colz");

  Beam2MoveX_bunch872Add->SetOption("colz");
  Beam2MoveY_bunch872Add->SetOption("colz");
  Beam1MoveX_bunch872Add->SetOption("colz");
  Beam1MoveY_bunch872Add->SetOption("colz");

  Beam2MoveX_bunch1783Add->SetOption("colz");
  Beam2MoveY_bunch1783Add->SetOption("colz");
  Beam1MoveX_bunch1783Add->SetOption("colz");
  Beam1MoveY_bunch1783Add->SetOption("colz");

  Beam2MoveX_bunch2063Add->SetOption("colz");
  Beam2MoveY_bunch2063Add->SetOption("colz");
  Beam1MoveX_bunch2063Add->SetOption("colz");
  Beam1MoveY_bunch2063Add->SetOption("colz");

  outFile->Write();



}
