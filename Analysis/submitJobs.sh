#!/bin/bash
workDir=`pwd`
mail="joscha.knolle@desy.de"
for j in configFill4954std.json #configFill4954prompt.json configFill4954loose.json configFill4954tight.json configFill4954verytight.json
#for j in configFill4937std.json
#for j in configFill4266std.json
do
    for m in SingleGaussUncorrelated #DoubleGauss SuperGauss
    do
        for c in 41 #281 872 1783 2063
        #for c in 81 875 1610 1690 1730
        #for c in 51 771 1631 2211 2674
        do
            for r in default #low high half onehalf
            do
                output="$workDir/log_${j%.*}_${m}_${c}_${r}.txt"
                qsub -j y -o $output -m eas -M $mail -l h_rt="5:00:00" jobShapeFitter.sh $workDir $j $m $c $r
            done
        done
    done
done
