#!/bin/bash
cd $1
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo "shapeFitter.py $2 -m $3 -c $4 -r $5"
python shapeFitter.py $2 -m $3 -c $4 -r $5
