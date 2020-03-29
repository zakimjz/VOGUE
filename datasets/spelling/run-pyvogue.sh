#!/bin/bash
vogue="../../vogue-model-wdur.py"
viterbi="../../vogue-wgap.py"
vghdir="$PWD"

FAM=("bad" "good") 
V=(10 10)
U=(5 5)

for i in 0 1 
do
    if [ ! -d $vghdir/${FAM[$i]} ]
    then
        mkdir $vghdir/${FAM[$i]}
    fi
    vghfile="$vghdir/${FAM[$i]}/${FAM[$i]}_v${V[$i]}_u${U[$i]}.vgh"
    train="$vghdir/${FAM[$i]}.train"
    alpha="$vghdir/SPELL.Alphabet"
    echo "$vogue -A $alpha -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile"
    $vogue -A $alpha -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile 
    for j in 0 1
    do
    	spell="$vghdir/${FAM[$j]}.test"
	outfile="$vghdir/${FAM[$i]}/${FAM[$j]}.out"
	echo "$viterbi -A $alpha -v $vghfile -t $spell > $outfile"
	$viterbi -A $alpha -v $vghfile -t $spell > $outfile
    done
done
