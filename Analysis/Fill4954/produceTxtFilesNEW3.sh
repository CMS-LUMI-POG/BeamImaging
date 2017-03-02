#!/bin/bash

rm -rf zeroBias*.txt
for i in {0..9}
do
##zero3
    for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias3/PCC_VdM_ZeroBias3_274100_ProMay312016_Event_AlwaysTrue/160531_211350/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias3/PCC_VdM_ZeroBias3_274100_ProMay312016_Event_AlwaysTrue/160531_211350/0000/$line" >> zeroBias3_$i.txt
    done
##zero2
    for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias2/PCC_VdM_ZeroBias2_274100_ProMay312016_Event_AlwaysTrue/160531_211331/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias2/PCC_VdM_ZeroBias2_274100_ProMay312016_Event_AlwaysTrue/160531_211331/0000/$line" >> zeroBias2_$i.txt
    done
##zero1
    for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias1/PCC_VdM_ZeroBias1_274100_ProMay312016_Event_AlwaysTrue/160531_211313/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias1/PCC_VdM_ZeroBias1_274100_ProMay312016_Event_AlwaysTrue/160531_211313/0000/$line" >> zeroBias1_$i.txt
    done
##zero4
    for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias4/PCC_VdM_ZeroBias4_274100_ProMay312016_Event_AlwaysTrue/160531_211413/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias4/PCC_VdM_ZeroBias4_274100_ProMay312016_Event_AlwaysTrue/160531_211413/0000/$line" >> zeroBias4_$i.txt
    done
##zero5
for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias5/PCC_VdM_ZeroBias5_274100_ProMay312016_Event_AlwaysTrue/160531_211430/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias5/PCC_VdM_ZeroBias5_274100_ProMay312016_Event_AlwaysTrue/160531_211430/0000/$line" >> zeroBias5_$i.txt
    done
##zero6
for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias6/PCC_VdM_ZeroBias6_274100_ProMay312016_Event_AlwaysTrue/160531_211447/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias6/PCC_VdM_ZeroBias6_274100_ProMay312016_Event_AlwaysTrue/160531_211447/0000/$line" >> zeroBias6_$i.txt
    done
##zero7
for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias7/PCC_VdM_ZeroBias7_274100_ProMay312016_Event_AlwaysTrue/160531_211507/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias7/PCC_VdM_ZeroBias7_274100_ProMay312016_Event_AlwaysTrue/160531_211507/0000/$line" >> zeroBias7_$i.txt
    done
##zero8
for line in `/afs/cern.ch/project/eos/installation/4.0.26-1/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias8/PCC_VdM_ZeroBias8_274100_ProMay312016_Event_AlwaysTrue/160531_211529/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/May272016_274100/ZeroBias8/PCC_VdM_ZeroBias8_274100_ProMay312016_Event_AlwaysTrue/160531_211529/0000/$line" >> zeroBias8_$i.txt
    done
done
