# VOGUE: Variable Order HMM with Duration

VOGUE is a variable order and gapped HMM with with duration. It uses sequence mining to extract frequent patterns in the data. It then uses these patterns to build a variable order HMM with explicit duration on the gap states, for sequence modeling and classification. VOGUE was applied to model protein sequences, as well as a number of other sequence datasets including weblogs.

Mohammed J. Zaki, Christopher D. Carothers, and Boleslaw K. Szymanski. VOGUE: a variable order hidden markov model with duration based on frequent sequence mining. ACM Transactions on Knowledge Discovery in Data, 4(1):Article 5, January 2010.


# How to

The datasets directory contains all the datasets used in the
paper, including the running example used in sec 2 and 3. Each
directory contains a shell script called run_pyvogue.sh that
contains the settings used to run VOGUE in the paper.

VOGUE accepts fasta format biological sequences or a spaced
format for other sequence data. See details below.

FASTA FORMAT SEQUENCES
----------------------

VOGUE accepts protein sequences in the fasta format as the
default. See the example training and test data in the prosite
and scop directory under datasets. 

With the default data (protein alphabet, fasta format)
Run VOGUE as follows:

    vogue-model-wdur.py -i TRAINFILE -s MINSUP -m MAXGAP -o VOGUEFILE

Run Viterbi as follows:

    vogue-wgap.py -v VOGUEFILE -t TESTFILE

See the example script in prosite and scop

SPACED FORMAT SEQUENCES
------------------------

* Sequence File:

VOGUE also accepts a fasta-like spaced format for modeling other
sequence data. You need to have a sequence data file in the
following format:

        > seq1
        xy abcz fgh abcefg
        > seq2 
        ghipr xy bbb fgh
        > .... 

Here > denotes a comment, followed by the sequence id or name,
and the next line gives the actual sequence consisting of space
separated words. Each word is treated as a symbol in the
alphabet. 

* Alphabet File:

In addition VOGUE requires an alphabet file, which lists all the
unique "symbols" in the alphabet. For the example above the file
would be as follows:

        xy abcz fgh abcefg ghipr bbb

That is the alphabet is given one a single line, separated by
spaces. See the example dataset used as the running example in
the paper (read the README file in the paper_example directory).
You can also look at the other datasets (cslogs, intrusion,
spelling)  

Now one can run VOGUE as follows:

    vogue-model-wdur.py -A ALPHABETFILE -i TRAINFILE -s MINSUP -m MAXGAP -o VOGUEFILE

Note the extra -A parameter that specifies the alphabet file.
There is no need to change the viterbi run. Run Viterbi as follows:

    vogue-wgap.py -v VOGUEFILE -t TESTFILE

See the example scripts in cslogs, intrusion, and spelling

PAPER RUNNING EXAMPLE:
-----------------------

The running example in the paper has only a single sequence,
therefore the use of -s 1 will find all sequences with
minimum support of 1. VOGUE allows one to prune based on
weighted support via the -w option. In this example, we set
-w 2 so that only those that occur 2 times across all
occurrences will be kept. In other words, regular support
(-s) increment the count only once per sequence, whereas
weighted support (-w) increments the count once per
occurrence.

OUTPUT FORMAT:
--------------

The viterbi output is as follows. For the running example it will
look like:

    test1 -8.53773680178 -21.610361858 -21.610361858 13.0726250562 13.0726250562 -8.53773680178 0.00249600410461 5

    The different fields are as follows:
    field 1: test sequence id (test1 in this case) 
    field 2: viterbi log probability
    field 3 and 4: null model log probability 
    field 5 and 6: log odd score (score-nullscore)
    field 7: viterbi score verified again
    field 8: time
    field 9: test sequence length

DEPENDENCIES:
--------------
VOGUE depends on the biopython package for fasta sequence parsing.
It also uses psyco to speed up the computation.



