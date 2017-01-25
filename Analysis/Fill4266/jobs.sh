#!/bin/bash

suffix=$1
outputDir=$3
workDirLocal=$2
workDirRemote=`pwd`

cd ${workDirLocal}/../
eval `scramv1 runtime -sh`
cd $workDirRemote

cp ${workDirLocal}/* .

root -l -b -q runOn_nTuples2.cc\(\"$suffix\"\)++
cmsStage *${suffix}.root $outputDir