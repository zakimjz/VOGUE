#!/usr/bin/env python3

import sys
from optparse import OptionParser
import vgs
import HmmerNull
import math
import time

# which alphabet to use, depends on order
AMINOALPHA = [HmmerNull.Amino, HmmerNull.Amino2, HmmerNull.Amino3,
              HmmerNull.Amino4]


class emission:
    def __init__(self):  # constructor
        self.B = []  # symbol esmission probabilities

    def emission_matrix(self, Sigma, itemcount, gapinfo, Qf, Qs,
                        maxgap, sym2idx, key2idx):

        self.B = [[]]*len(key2idx)
        # set the emission probs for first states
        for first in Qf:
            b = [0.0]*len(sym2idx)  # init to all zeros
            b[sym2idx[first]] = 1.0
            self.B[key2idx['first', first]] = b

        # set the emission probs for second states
        for second in Qs:
            b = [0.0]*len(sym2idx)  # init to all zeros
            b[sym2idx[second]] = 1.0
            self.B[key2idx['second', second]] = b

        frac = 0.01
        # set the emission probs for universal state
        b = [frac/len(Sigma)]*len(sym2idx)  # init to frac/|Sigma|
        bsum = 0
        for item in itemcount:
            bsum += itemcount[item]
        for item in itemcount:
            b[sym2idx[item]] += ((1-frac)*itemcount[item])/bsum
        self.B[key2idx['universal']] = b

        # set the emission probs for gap states
        for pair in gapinfo:
            bsum = 0.0
            # first compute sum of all gap symbols in the 2-seq given by pair
            for k in gapinfo[pair][2]:
                bsum += gapinfo[pair][2][k]
            b = [0.0]*len(sym2idx)  # output symbol hash

            if not gapinfo[pair][2]:  # dict is empty
                for s in Sigma:
                    b[sym2idx[s]] = 1.0/len(Sigma)
            else:
                for s in Sigma:
                    b[sym2idx[s]] = frac/len(Sigma)

            for k in gapinfo[pair][2]:
                b[sym2idx[k]] += ((1-frac)*gapinfo[pair][2][k])/bsum

            if options.nodur:
                # for all gaps in 1..maxgap set the output emission prob
                for g in range(1, maxgap+1):
                    self.B[key2idx['gap', pair, g]] = b
            else:
                # set the output emission prob
                self.B[key2idx['gap', pair]] = b

        return self.B


class transition:
    def __init__(self):
        self.pi = []  # init state probs
        self.A = []  # transition matrix A

    def transition_matrix(self, gapinfo, Qf, Qs, Qu, Ng, maxgap, key2idx):
        # init state probs
        # only to first states or universal gap
        self.pi = [0.0]*len(key2idx)
        for s in Qs:
            self.pi[key2idx['second', s]] = 0.0
        for key in gapinfo:
            if options.nodur:
                for g in range(1, maxgap+1):
                    self.pi[key2idx['gap', key, g]] = 0.0
            else:
                self.pi[key2idx['gap', key]] = 0.0
        frac = 0.01
        for f in Qf:
            self.pi[key2idx['first', f]] = ((1-frac)*Qf[f])/Qu
        self.pi[key2idx['universal']] = frac

        self.A = [[]]*len(key2idx)
        # transitions from first states
        # only to second states or first gap state per freq 2-seq
        for f in Qf:
            #print "processing first", f
            Af = [0.0]*len(key2idx)
            for f2 in Qf:
                Af[key2idx['first', f2]] = 0.0
            Af[key2idx['universal']] = 0.0
            for s in Qs:
                Af[key2idx['second', s]] = 0.0
            for key in gapinfo:
                if options.nodur:
                    for g in range(1, maxgap+1):
                        Af[key2idx['gap', key, g]] = 0.0
                else:
                    Af[key2idx['gap', key]] = 0.0
                if key[0] == f:  # first state matches key
                    gsum = sum(gapinfo[key][1][1:])  # skip gap=0
                    Af[key2idx['second', key[1]]
                       ] = gapinfo[key][1][0]/Qf[f]
                    if (maxgap > 0):
                        if options.nodur:
                            Af[key2idx['gap', key, 1]] = gsum/Qf[f]
                        else:
                            Af[key2idx['gap', key]] = gsum/Qf[f]
            self.A[key2idx['first', f]] = Af

        # transitions from second states are same as pi
        for s in Qs:
            self.A[key2idx['second', s]] = self.pi

        # transitions from universal gap state
        # only to self or first states
        frac = 0.01
        Au = [0.0]*len(key2idx)
        #Au[key2idx['universal']] = 0.01;
        Au[key2idx['universal']] = frac
        for s in Qs:
            Au[key2idx['second', s]] = 0.0
        for f in Qf:
            Au[key2idx['first', f]] = ((1-frac)*Qf[f])/Qu
        for key in gapinfo:
            if options.nodur:
                if (maxgap > 0):
                    Au[key2idx['gap', key, 1]] = frac/Ng
                for g in range(2, maxgap+1):
                    Au[key2idx['gap', key, g]] = 0.0
            else:
                #Au[key2idx['gap',key]] = (frac-0.01)/Ng
                Au[key2idx['gap', key]] = 0.0

        self.A[key2idx['universal']] = Au

        # transitions from gap states
        # only to next gap state or to second state
        if options.nodur:
            for key in sorted(gapinfo):
                for g in range(1, maxgap+1):
                    #print "processing gap", key, g
                    Ag = [0.0]*len(key2idx)
                    # set all to zero
                    for f in Qf:
                        Ag[key2idx['first', f]] = 0.0
                    Ag[key2idx['universal']] = 0.0
                    for key2 in gapinfo:
                        for g2 in range(1, maxgap+1):
                            Ag[key2idx['gap', key2, g2]] = 0.0
                    for s in Qs:
                        Ag[key2idx['second', s]] = 0.0

                    # now just modify the trans to next gap or second state
                    if g == maxgap:
                        Ag[key2idx['second', key[1]]] = 1.0
                    else:
                        gsum = sum(gapinfo[key][1][1:])
                        fsum = sum(gapinfo[key][1][g+1:])
                        tprob = 0.0
                        if gsum > 0:
                            tprob = fsum/gsum
                        else:
                            tprob = 0.0
                        Ag[key2idx['gap', key, g+1]] = tprob
                        Ag[key2idx['second', key[1]]] = 1 - tprob
                    self.A[key2idx['gap', key, g]] = Ag
        else:
            for key in sorted(gapinfo):
                #print "processing gap", key, g
                Ag = [0.0]*len(key2idx)
                # set all to zero
                for f in Qf:
                    Ag[key2idx['first', f]] = 0.0
                Ag[key2idx['universal']] = 0.0
                for key2 in gapinfo:
                    Ag[key2idx['gap', key2]] = 0.0
                for s in Qs:
                    Ag[key2idx['second', s]] = 0.0
                # now just modify the trans to second state
                Ag[key2idx['second', key[1]]] = 1.0

                self.A[key2idx['gap', key]] = Ag


class duration:
    def __init__(self):
        self.rho = []  # state dur probs

    def duration_matrix(self, gapinfo, Qf, Qs, Qu, Ng, maxgap, key2idx):
        self.rho = [[]]*len(key2idx)
        for f in Qf:
            self.rho[key2idx['first', f]] = [0.0]*(maxgap+1)
            self.rho[key2idx['first', f]][0] = 1.0

        self.rho[key2idx['universal']] = [0.0]*(maxgap+1)
        self.rho[key2idx['universal']][0] = 1.0
        for s in Qs:
            self.rho[key2idx['second', s]] = [0.0]*(maxgap+1)
            self.rho[key2idx['second', s]][0] = 1.0

        for key in sorted(gapinfo):
            self.rho[key2idx['gap', key]] = [0.0]*(maxgap+1)
            self.rho[key2idx['gap', key]][0] = 1.0
            for g in range(1, maxgap+1):
                tsum = sum(gapinfo[key][1][1:])
                gsum = gapinfo[key][1][g]
                tprob = 0.0
                if tsum > 0:
                    tprob = gsum*1.0/tsum
                self.rho[key2idx['gap', key]][g-1] = tprob


def check_sum_eq_one(T):
    if type(T) == dict:
        sary = [T[x] for x in T]
    else:
        sary = T
    tsum = sum(sary)
    if math.fabs(1.0-tsum) < 10e-5:
        return True
    else:
        print("SUM: ", tsum)
        if type(T) == dict:
            for x in T:
                T[x] = T[x]/tsum
        else:
            for i in range(len(T)):
                T[i] = T[i]/tsum
        return True


def read_alphabet(afile):
    S = []
    f = open(afile, "r")
    # assume space separated words
    for line in f:
        for word in line.strip().split():
            S.append(word)
    f.close()
    return S


class MyOptionParser (OptionParser):
    def check_required(self, *opt):
        for op in opt:
            option = self.get_option(op)
            # Assumes the option's 'default' is set to None!
            if getattr(self.values, option.dest) is None:
                self.error("%s option not supplied" % option)


parser = MyOptionParser()
# help message is automatically provided
parser.add_option('-A', '--AfromFile', dest='AfromFile',
                  help='read alphabet from given file')
parser.add_option('-a', '--alphabet', dest='alphabet',
                  help='which alphabet to use (protein, kvogue, dna)')
parser.add_option('-c', '--artifacts', type='float',
                  dest='artifacts', help='remove artifacts; thresh ratio')
parser.add_option('-d', '--nodur', action='store_true',
                  dest='nodur', help='no durations')
parser.add_option('-i', '--infile', dest='infile',
                  help='input file')
parser.add_option('-f', '--format', dest='format',
                  help='file format string (e.g., fasta)')
parser.add_option('-m', '--maxgap', type='int', dest='maxgap',
                  help='maximum gap between two symbols')
parser.add_option('-n', '--trainsize', type='int', dest='trainsize',
                  help='number of training sequences to use')
parser.add_option('-o', '--output', dest='outfile',
                  help='print VOGUE to a file')
parser.add_option('-R', '--hmmorder', type='int', dest='hmmorder',
                  help='order of HMM (#chars/symbol)')
parser.add_option('-s', '--minsup', type='int', dest='minsup',
                  help='absolute minsup value (#of sequences)')
parser.add_option('-S', '--minsupratio', type='float', dest='minsupratio',
                  help='relative minsup value (%of sequences)')
parser.add_option('-w', '--wminsup', type='int', dest='wminsup',
                  help='weighted minsup value (#of occurrences)')
parser.add_option('-v', '--verify', action='store_true', dest='verify',
                  help='verify if probs add to 1.0')
options, args = parser.parse_args(sys.argv[1:])
parser.check_required('-i', '-m')

t0 = time.time()

if options.verify is None:
    options.verify = True

if options.hmmorder is None:
    options.hmmorder = 1  # default one char/symbol

if options.format is None:
    options.format = "fasta"

if options.alphabet is None:
    options.alphabet = "protein"

kvogue = False
if options.alphabet == "protein":
    Sigma = AMINOALPHA[options.hmmorder-1]
elif options.alphabet == "kvogue":
    Sigma = HmmerNull.KVOGUE
    kvogue = True
elif options.alphabet == "dna":
    Sigma = HmmerNull.DNA

if options.AfromFile is not None:
    Sigma = read_alphabet(options.AfromFile)
    options.format = "spaced"

if options.outfile is not None:
    outf = open(options.outfile, 'w')
else:
    outf = sys.stdout

#print options.infile, options.maxgap, options.minsup,\
#    options.outfile, options.wminsup


itemcount, gapinfo = vgs.vgs(options.infile, options.maxgap,
                             options.minsup, options.minsupratio,
                             options.wminsup, options.format,
                             kvogue, options.hmmorder,
                             options.artifacts,
                             options.trainsize)

# pprint(itemcount)
# pprint(gapinfo)

Qu = 0.0  # universal gap, sum of all 2-freq weighted sups
Qf = {}  # first states, keep freq info
Qs = {}  # second states
Ng = len(gapinfo)  # number of keys = # of 2-seq

for key in gapinfo:
    if key[0] not in Qf:
        Qf[key[0]] = 0.0
    Qf[key[0]] += gapinfo[key][0]  # total freq at first state
    Qs[key[1]] = True  # which are the second states
    Qu += gapinfo[key][0]  # total freq at universal gap

#print "VOGUE States: ", len(Qf), len(Qs), Ng

# create a symbol to index map, to use list instead of dict
sym2idx = {}
idx = 0
for s in Sigma:
    sym2idx[s] = idx
    idx += 1

# create a key to index map, so we can use list instead of dict
# for both emissions and transition probs
key2idx = {}  # hash map from key to idx
idx = 0
for f in Qf:
    key2idx['first', f] = idx
    idx += 1

for key in sorted(gapinfo):
    if options.nodur:
        for g in range(1, options.maxgap+1):
            key2idx['gap', key, g] = idx
            idx += 1
    else:
        key2idx['gap', key] = idx
        idx += 1

for s in Qs:
    key2idx['second', s] = idx
    idx += 1
key2idx['universal'] = idx
#print "TOTAL VOGUE STATES ", len(key2idx)

E = emission()
E.emission_matrix(Sigma, itemcount, gapinfo, Qf, Qs, options.maxgap,
                  sym2idx, key2idx)

T = transition()
T.transition_matrix(gapinfo, Qf, Qs, Qu, Ng, options.maxgap, key2idx)

if options.nodur is None:
    D = duration()
    D.duration_matrix(gapinfo, Qf, Qs, Qu, Ng, options.maxgap, key2idx)

if options.verify:
    if not check_sum_eq_one(T.pi):
        print("PROBLEM: T.pi = %s" % (T.pi))
        sys.exit(-1)
    for key in E.B:
        if not check_sum_eq_one(key):
            print("PROBLEM: E.B= % s" % (key))
            sys.exit(-1)
    for i, key in enumerate(T.A):
        if not check_sum_eq_one(key):
            print("PROBLEM: T.A %d = %s" % (i, key))
            sys.exit(-1)

if options.nodur:
    print("CompactModel", file=outf)
else:
    print("DurationCompactModel", file=outf)
print("M:", len(sym2idx), file=outf)
print("N:", len(key2idx), file=outf)
print("N_1:", len(Qf), file=outf)

if options.nodur:
    print("N_2:", len(key2idx) - len(Qf) - len(Qs), file=outf)
else:
    print("N_2:", len(key2idx) - len(Qf) - len(Qs) - 1, file=outf)
print("N_3:", len(Qs), file=outf)

if options.nodur:
    if options.maxgap == 0:
        print("MAXGAP:", options.maxgap+1, file=outf)
    else:
        print("MAXGAP:", options.maxgap, file=outf)
else:
    print("MAXGAP: 1", file=outf)

print("A:", file=outf)
for q, state_probs in enumerate(T.A):
    for i, v in enumerate(state_probs):
        if v > 0.0:
            print("%d %d %g" % (q, i, v), file=outf)

print("B:", file=outf)
for q, emit_probs in enumerate(E.B):
    for i, v in enumerate(emit_probs):
        if v > 0.0:
            print("%d %d %g" % (q, i, v), file=outf)

if options.nodur is None:
    print("rho:", file=outf)
    for q, dur_probs in enumerate(D.rho):
        for i, p in enumerate(dur_probs):
            if p > 0.0:
                print("%d %d %g" % (q, i, p), file=outf)

print("pi:", file=outf)
for i, p in enumerate(T.pi):
    if p > 0.0:
        print("%d %g" % (i, p), file=outf)

t1 = time.time()
print("total time", t1-t0)
