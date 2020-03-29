#!/bin/bash
vogue="../../vogue-model-wdur.py"
viterbi="../../vogue-wgap.py"
vghdir="$PWD"

FAM=("USER0" "USER1" "USER2" "USER3" "USER4" "USER5" "USER6" "USER7"
    "USER8")
V=(2 2 2 2 2 2 2 2 2)
U=(1 1 1 1 1 1 1 1 1)

for i in 0 1 2 3 4 5 6 7 8 
do
    if [ ! -d $vghdir/${FAM[$i]} ]
    then
        mkdir $vghdir/${FAM[$i]}
    fi
    vghfile="$vghdir/${FAM[$i]}/${FAM[$i]}_v${V[$i]}_u${U[$i]}.vgh"
    train="$vghdir/${FAM[$i]}.train"
    alpha="$vghdir/USERS.Alphabet"
    echo "$vogue -A $alpha -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile"
    $vogue -A $alpha -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile 
    for j in 0 1 2 3 4 5 6 7 8
    do
    	intrusion="$vghdir/${FAM[$j]}.test"
	outfile="$vghdir/${FAM[$i]}/${FAM[$j]}.out"
	echo "$viterbi -A $alpha -v $vghfile -t $intrusion > $outfile"
	$viterbi -A $alpha -v $vghfile -t $intrusion > $outfile
    done
done
