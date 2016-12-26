#!/bin/bash

#outputDir=/store/user/jsalfeld/vdmAnalysis/beamImaging/statisticsEffect/
#outputDir=/store/user/jsalfeld/vdmAnalysis/beamImaging/statisticsEffect_loosenedConstraints2/

#outputDir=/store/user/jsalfeld/vdmAnalysis/beamImaging/noCorrelations2/

outputDir=/store/user/jsalfeld/vdmAnalysis/beamImaging/noCorrelations2_noOverlapInWidth/

workDirLocal=`pwd`
echo $workDirLocal

script=jobs.sh

for i in {2000..3000}
do
bsub -o out.%J -q 2nw  $script $i $workDirLocal $outputDir
done