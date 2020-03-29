#!/bin/bash
vogue="../../vogue-model-wdur.py"
viterbi="../../vogue-wgap.py"
vghdir="$PWD"

FAM=("F1-49417" "F2-46458" "F3-46626" "F4-46689" "F5-46997"
	"F6-47095" "F7-47113" "F8-48508" "F9-69118" "F10-81296") 
V=(1 1 1 1 1 1 1 1 1 1)
U=(30 30 30 30 30 30 30 30 30 30)

for i in 0 1 2 3 4 5 6 7 8 9
do
    if [ ! -d $vghdir/${FAM[$i]} ]
    then
        mkdir $vghdir/${FAM[$i]}
    fi
    vghfile="$vghdir/${FAM[$i]}/${FAM[$i]}_v${V[$i]}_u${U[$i]}.vgh"
    train="$vghdir/${FAM[$i]}.train"
    echo "$vogue -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile"
    $vogue -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile 
    for j in 0 1 2 3 4 5 6 7 8 9
    do
    	scop="$vghdir/${FAM[$j]}.test"
	outfile="$vghdir/${FAM[$i]}/${FAM[$j]}.out"
	echo "$viterbi -v $vghfile -t $scop > $outfile"
	$viterbi -v $vghfile -t $scop > $outfile
    done
done
