#!/bin/bash
vogue="../../vogue-model-wdur.py"
viterbi="../../vogue-wgap.py"
vghdir="$PWD"

FAM=(00662 00670 00561 00064 00154 00224 00271 00343 00397 00443)
V=(1 1 1 1 1 1 1 1 1 1)
U=(1 1 1 1 2 2 8 1 1 1)

for i in 0 1 2 3 4 5 6 7 8 9 
do
    if [ ! -d $vghdir/PDOC${FAM[$i]} ]
    then
        mkdir $vghdir/PDOC${FAM[$i]}
    fi
    vghfile="$vghdir/PDOC${FAM[$i]}/PDOC${FAM[$i]}_v${V[$i]}_u${U[$i]}.vgh"
    train="$vghdir/class_PDOC${FAM[$i]}_positive.train"
    echo "$vogue -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile"
    $vogue -i $train -s ${V[$i]} -m ${U[$i]} -o $vghfile 
    for j in 0 1 2 3 4 5 6 7 8 9
    do
    	prosite="$vghdir/class_PDOC${FAM[$j]}_positive.test"
	outfile="$vghdir/PDOC${FAM[$i]}/PDOC${FAM[$j]}.out"
	echo "$viterbi -v $vghfile -t $prosite > $outfile"
	$viterbi -v $vghfile -t $prosite > $outfile
    done
done
