# BeamImaging

This is the repository for code related to 2 dimensional proton bunch densitiy measurements.

* Setting up the CMSSW environment should be done before checking out. It has been tested to work in CMSSW_7_4_2.

* The code is organized in terms of fills. "Analysis/FillWXYZ" contains a scripts "runOn_nTuples2.cc" to produce 2d histograms, and the scripts "jobs.sh" and "submitVdMjobs.sh" are helpful to run over the nTuples and batch jobs submission.

* The directory "beamdensities" contains customized RooFit PDF classes which have to be compiled using "makeFile3". The PDF classes are then used by the scripts "Analysis/FillWXYZ/vdmScanTreeAnalyzerABC_onData2_studyUnc.C" and have to be loaded before they can be used:

```C
    root -l
    gSystem->Load("../../beamdensities/MyPdfXYZ.so")
    .x vdmScanTreeAnalyzerDG_onData2_studyUnc.C++
```

* The fits of the vertex distributions to the beam shape models are done with the scripts "Analysis/FillWXYZ/vdmScanTreeAnalyzerABC_onData2_studyUnc.C" where each fit model has a separate file. Some of these files were changed to run over one bunch crossing specified as argument only, while the others run without any argument.

* VdM scan measurements of the overlap integral simulations and correction factor calculations are done with the scripts "overlapDiff_AB.C" where each class of fit models has a separate file. The function takes as argument a string of the bunch crossing and the fit model concatenated.

* Some python scripts just extract values from ROOT files and output them to the screen in CSV style: "computeChiSquares.py" (chi2 values for all models and bcids), "gatherFromToys.py" (correction factors for all models and bcids), "gatherParameters" (Gaussian parameters for all models and bcids).

* Some python scripts just make PAS-style plots of the results: "summaryPlots.py" (all relevant Beam Imaging plots), "plotOtherCoordinate.py" (plot of x mean and width during y scan and vice versa)
