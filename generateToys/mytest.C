void
mytest(TString const& _input)
{
  gSystem->Load("MyPdfV4.so");
  gSystem->Load("MyPdfV3.so");
  gROOT->LoadMacro("vdmScanTreeAnalyzerDG_toys.C+");

  gROOT->ProcessLine("vdmScanTreeAnalyzerDG_toys(\"" + _input + "\")");
}
