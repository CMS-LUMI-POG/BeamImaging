#!/bin/bash

#outputDir=/store/user/jsalfeld/vdmAnalysis/beamImaging/statisticsEffect/
#outputDir=/afs/cern.ch/work/j/jsalfeld/MIT/vdmAnalysis/CMSSW_7_4_2/src/vdm_2582015_Analysis/outFiles2//
#outputDir=/store/user/jsalfeld/vdmScan_2482025/BI_histosNEW7_ext/
#outputDir=/store/user/jsalfeld/vdmScan_2482025/BI_histosNEW7_ext_withLSC_z/

#outputDir=/store/user/jsalfeld/vdmScanPromptReco2015NEW
outputDir=/store/group/comm_luminosity/jsalfeld/BI_PromptReco2015/
workDirLocal=`pwd`
echo $workDirLocal

script=jobs.sh

for i in {1..8}
do
    for j in {0..9}
    do
	bsub -o out.%J -q cmscaf1nd  $script zeroBias${i}_${j}.txt $workDirLocal $outputDir
    done
done