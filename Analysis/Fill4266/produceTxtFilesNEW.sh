#!/bin/bash

rm -rf zeroBias*.txt
for i in {0..9}
do
##zero3
    for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias3/PCC_ZeroBias3_VdMScans_150825_nVtxFix/150831_134550/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias3/PCC_ZeroBias3_VdMScans_150825_nVtxFix/150831_134550/0000/$line" >> zeroBias3_$i.txt
    done
##zero2
    for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias2/PCC_ZeroBias2_VdMScans_150825_nVtxFix/150831_134510/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias2/PCC_ZeroBias2_VdMScans_150825_nVtxFix/150831_134510/0000/$line" >> zeroBias2_$i.txt
    done
##zero1
    for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias1/PCC_ZeroBias1_VdMScans_150825_nVtxFix/150831_134443/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias1/PCC_ZeroBias1_VdMScans_150825_nVtxFix/150831_134443/0000/$line" >> zeroBias1_$i.txt
    done
##zero4
for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias4/PCC_ZeroBias4_VdMScans_150825_nVtxFix/150831_134629/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias4/PCC_ZeroBias4_VdMScans_150825_nVtxFix/150831_134629/0000/$line" >> zeroBias4_$i.txt
    done
##zero5
for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias5/PCC_ZeroBias5_VdMScans_150825_nVtxFix/150831_134649/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias5/PCC_ZeroBias5_VdMScans_150825_nVtxFix/150831_134649/0000/$line" >> zeroBias5_$i.txt
    done
##zero6
for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias6/PCC_ZeroBias6_VdMScans_150825_nVtxFix/150831_134708/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias6/PCC_ZeroBias6_VdMScans_150825_nVtxFix/150831_134708/0000/$line" >> zeroBias6_$i.txt
    done
##zero7
for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias7/PCC_ZeroBias7_VdMScans_150825_nVtxFix/150831_134731/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias7/PCC_ZeroBias7_VdMScans_150825_nVtxFix/150831_134731/0000/$line" >> zeroBias7_$i.txt
    done
##zero8
for line in `/afs/cern.ch/project/eos/installation/0.3.35/bin/eos.select ls /store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias8/PCC_ZeroBias8_VdMScans_150825_nVtxFix/150831_134748/0000/*$i.root`
    do
	echo "root://eoscms//eos/cms/store/group/comm_luminosity/PCC/VdM/0150825/ZeroBias8/PCC_ZeroBias8_VdMScans_150825_nVtxFix/150831_134748/0000/$line" >> zeroBias8_$i.txt
    done
done