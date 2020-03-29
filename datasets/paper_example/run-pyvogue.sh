#!/bin/bash
vogue="../../vogue-model-wdur.py"
viterbi="../../vogue-wgap.py"
vghdir="$PWD"

vghfile="$vghdir/s1w2m2.vgh"
train="$vghdir/train.fasta"
alpha="$vghdir/train.alpha"
echo "$vogue -i $train -s 1 -w 2  -m 2 -o $vghfile"
$vogue -i $train -A $alpha -s 1 -w 2 -m 2 -o $vghfile 
    	
test="$vghdir/test.fasta"
echo "$viterbi -v $vghfile -t $test"
$viterbi -v $vghfile -t $test
