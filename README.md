# BeamImaging

This is the repository for code related to 2 dimensional proton bunch densitiy measurements.

Setting up the CMSSW environment should be done before checking out. It has been tested to work in CMSSW_7_4_2.

The code is organized in terms of fills. "Analysis/FillWXYZ" contains a scripts "runOn_nTuples2.cc" to produce 2d histograms, and the scripts "jobs.sh" and "submitVdMjobs.sh" are helpful to run over the nTuples and batch jobs submission.

The directory "beamdensities" contains customized RooFit PDF classes which have to be compiled using "makeFile3".

The PDF classes are then used by the scripts "Analysis/FillWXYZ/vdmScanTreeAnalyzerDG_onData2_studyUnc.C" to fit the 2d histograms. The library has to tbe loaded before it can be used:

root -l
gSystem->Load("../../beamdensities/MyPdfXYZ.so")
.x vdmScanTreeAnalyzerDG_onData2_studyUnc.C++

There are other scripts to be added for the future for simulation and proper treatment of the bunch profiles.