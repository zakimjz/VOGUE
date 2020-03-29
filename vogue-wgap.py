#!/usr/bin/env python
import sys, optparse, math
from optparse import OptionParser
import HmmerNull
import time
import ParseSeq
import psyco
from psyco.classes import *

verify_score = True

#which alphabet to use, depends on order
AMINOALPHA = [HmmerNull.Amino, HmmerNull.Amino2, HmmerNull.Amino3,
         HmmerNull.Amino4]

class HMM():
    def __init__(self, voguefile, Sigma):
        self.A = {} #transition probs, init probs as A[0,x]
        self.B = {} #emission probs
        self.D = {} #duration probs
        self.m = self.n = self.n1 = self.n2 = self.n3 = self.maxgap = 0
        self.read_vgh_file(voguefile, Sigma)
    
    def read_vgh_file(self,voguefile, Sigma):
        f = open(voguefile, "r")
        firstline = f.readline().strip()
        compactmodel = False
        durationcompactmodel = False
        if firstline == "CompactModel":
            compactmodel = True
        elif firstline == "DurationCompactModel":
            durationcompactmodel = True
        
        if compactmodel or durationcompactmodel:
            self.m = int(f.readline().strip().split()[1])
        else:
            self.m = int(firstline.split()[1])
        
        self.n = int(f.readline().strip().split()[1])
        self.n1 = int(f.readline().strip().split()[1])
        self.n2 = int(f.readline().strip().split()[1])
        self.n3 = int(f.readline().strip().split()[1])
        self.maxgap = int(f.readline().strip().split()[1])
        
        if compactmodel or durationcompactmodel:
            self._read_compact_ABPD(f, Sigma, durationcompactmodel)
        else:
            self._read_ABPD(f, Sigma)

        f.close()
        self._transitions_to = self._calculate_to_transitions()
        self.n +=1 # increase state to include beg state

    def _read_compact_ABPD(self, f, Sigma, duration):
        s = f.readline().strip() #this is A:
        if (s != "A:"): sys.exit(-1)

        inline = f.readline().strip().split()
        while(inline[0] != "B:"):
            i = int(inline[0])
            j = int(inline[1])
            p = float(inline[2])
            if p > 0.0:
                self.A[i+1,j+1] = math.log(p,2)
            inline = f.readline().strip().split()
        
        if (inline[0] != "B:"): sys.exit(-1)
        inline = f.readline().strip().split()
        if duration:
            sym = "rho:"
        else: sym = "pi:"
        while inline[0] != sym:
            i = int(inline[0])
            j = int(inline[1])
            p = float(inline[2])
            if p > 0.0:
                self.B[i+1,Sigma[j]] = math.log(p,2)
            inline = f.readline().strip().split()
        
        if duration:
            if (inline[0] != "rho:"): sys.exit(-1)
            inline = f.readline().strip().split()
            while inline[0] != "pi:":
                i = int(inline[0])
                j = int(inline[1])
                p = float(inline[2])
                self.D[i+1,j+1] = math.log(p,2)
                inline = f.readline().strip().split()

            self.D[0,1] = 0.0 #init state has duration 1
        else:
            #each state has duration 1
            self.D[0,1] = 0.0 
            for i in range(self.n):
                self.D[i+1,1] = 0.0 

        if (inline[0] != "pi:"): sys.exit(-1)
        inline = f.readline().strip().split()
        while(len(inline) > 0):
            i = int(inline[0])
            p = float(inline[1])
            if p > 0.0:
                self.A[0,i+1] = math.log(p,2)
            inline = f.readline().strip().split()


    def _read_ABPD(self, f, Sigma):
        s = f.readline().strip() #this is A:
        if (s != "A:"): sys.exit(-1)
        for i in range(self.n):
            probs = [float(p) for p in f.readline().strip().split()]
            for j, p in enumerate(probs):
                if p > 0.0:
                    self.A[i+1,j+1] = math.log(p,2)
        
        s = f.readline().strip() #this is B:
        if (s != "B:"): sys.exit(-1)
        
        for i in range(self.n):
            emits = [float(p) for p in f.readline().strip().split()]
            for j, p in enumerate(emits):
                if p > 0.0:
                    self.B[i+1,Sigma[j]] = math.log(p,2)
        
        s = f.readline().strip() #this is D:
        if (s != "rho:"): sys.exit(-1)
        
        #these are the durations for each state.
        for i in range(self.n):
            rhos = [float(p) for p in f.readline().strip().split()]
            for j, p in enumerate(rhos):
                if p > 0.0:
                    self.D[i+1,j+1] = math.log(p,2)
        self.D[0,1] = 0.0 #init state has duration 1

        s = f.readline().strip() #this is pi:
        if (s != "pi:"): sys.exit(-1)
        
        pi = [float(p) for p in f.readline().strip().split()]
        for i, p in enumerate(pi):
                if p > 0.0:
                    self.A[0,i+1] = math.log(p,2)



    
    def _calculate_to_transitions(self):
        to_transitions = {}
        
        # loop over all of the different transitions
        for trans_key in self.A.keys():
            # if the letter to 'transition from' already exists, add the
            # new letter which can be 'transitioned to' to the list
            try:
                to_transitions[trans_key[1]].append(trans_key[0])
            # otherwise create the list and add the letter
            except KeyError:
                to_transitions[trans_key[1]] = []
                to_transitions[trans_key[1]].append(trans_key[0])
       
        #sort the lists
        for k in to_transitions:
            to_transitions[k].sort()

        return to_transitions
    
    def transitions_to(self, state):
        """Get all transitions which can happen from the given state.
        
        This returns all letters which the given state_letter is allowed
        to transition to. An empty list is returned if no letters are possible.
        """
        try:
            return self._transitions_to[state]
        except KeyError:
            return []

    
    def _trace_back(self, viterbi_probs, pred_state_seq, seqpos):
        #print seqpos
        #print "VP:"
        #print_dict(viterbi_probs)
        #print "PSS:"
        #print_dict(pred_state_seq)
        
        # find the max prob at last pos
        all_probs = {}
        for state in range(self.n):
            if (state,seqpos) in viterbi_probs:
                all_probs[state] = viterbi_probs[(state,seqpos)]
        
        if len(all_probs) > 0:
            state_path_prob = max(all_probs.values())
        
        # find the last pointer we need to trace back from
        last_state = None
        for state in all_probs.keys():
            if all_probs[state] == state_path_prob:
                last_state = state
        
        assert last_state is not None, "Can't trace from last state!"
        
        #print "MAX PROB", state_path_prob, last_state, seqpos
        
        # --- traceback
        traceback_seq = []
        
        cur_state = last_state
        pos = seqpos
        while (cur_state, pos) in pred_state_seq:
            prev_state, gap = pred_state_seq[(cur_state, pos)]
            for p in range(pos-gap, pos):
                traceback_seq.append(cur_state)
            pos = pos-gap
            cur_state = prev_state

        assert cur_state == 0, "traceback failed" + str(cur_state) +\
                " " + str(pos)
        traceback_seq.append(cur_state)
        
        # put the traceback sequence in the proper orientation
        traceback_seq.reverse()
        
        return traceback_seq, state_path_prob
    
    def viterbi(self, sequence):
        """Calculate the most probable state path using the Viterbi algorithm.
        Arguments:
        
        o sequence -- A Seq object with the emission sequence that we
        want to decode.
        """
        
        viterbi_probs = {}
        pred_state_seq = {}
        viterbi_probs[0, -1] = 0.0 #initial state is 0
        
        # --- recursion
        # loop over the training sequence (i = 1 .. L)
        for i in range(len(sequence)):
            #print "POS: ", i, sequence[i]
            # now loop over all of the letters in the state path
            somestatefound = False
            for cur_state in range(self.n):
                emission_part = {}
                for d in range(1,self.maxgap+1):
                    if (cur_state,d) in self.D:
                        emitprob = self.get_subseq_prob(sequence, cur_state, i, d)
                        if emitprob is not None: 
                            emission_part[d] = emitprob
                        else:
                            break

                if len(emission_part) == 0:
                    continue #skip to next state
                #print "CURRENT STATE:", cur_state
                
                # loop over all possible states
                prob_ary = {}
                for prev_state in self.transitions_to(cur_state):
                    for d in emission_part:
                        if (prev_state, i-d) in viterbi_probs:
                            #print "PREV STATE:", prev_state, d
                            trans_part = self.A[prev_state, cur_state]
                            emit_part = emission_part[d] 
                            viterbi_part = viterbi_probs[prev_state, i-d]
                            prev_prob = viterbi_part + trans_part + emit_part +\
                                self.D[cur_state,d]
                            prob_ary[prev_state,d] = prev_prob
                            #print "PA: ", "(", prev_state, d, ") ",\
                                    #            prev_prob, cur_state 
                
                
                # finally calculate the viterbi probability using the max
                if len(prob_ary) > 0:
                    somestatefound = True
                    max_prob, max_state, gap = self.get_max_arg(prob_ary)
                    viterbi_probs[(cur_state, i)] = max_prob
                    pred_state_seq[(cur_state, i)] = (max_state, gap)

                    #print "VP:", cur_state, i, viterbi_probs[(cur_state, i)] 
                    #print "PSS:", cur_state, i, pred_state_seq[(cur_state, i)] 
            
            assert somestatefound, "no state found for pos " + str(i)
        
        return self._trace_back(viterbi_probs, pred_state_seq,
                                len(sequence)-1)

    
    def get_subseq_prob(self, sequence, state, i, d):
        sum = 0.0
        for pos in range(i-d+1, i+1):
            if pos < len(sequence) and pos >= 0:
                if (state, sequence[pos]) in self.B:
                    sum += self.B[state, sequence[pos]]
                else:
                    return None
            else: return None
        return sum
    
    """ return the maximum value and the argmax """
    def get_max_arg(self,probs):
        first = True
        for (state,gap) in probs:
            if first:
                first = False
                mval = probs[state,gap]
                mstate = state
                mgap = gap

            if probs[state,gap] > mval:
                mval = probs[state,gap]
                mstate = state
                mgap = gap

        return mval, mstate, mgap

    
    def write_model(self):
        print self.m, self.n, self.n1, self.n2, self.n3, self.maxgap
        print "B:"
        for key in sorted(self.B):
            print "%s: %f" % (key, self.B[key])
        print "A:"
        for key in sorted(self.A):
            print "%s: %4.3f" % (key, self.A[key])
        print "D:"
        for key in sorted(self.D):
            print "%s: %4.3f" % (key, self.D[key])



def verify_score(v, A, B, D, seq, laststate):
    prob = 0.0
    prev = 0 #init state is 0 
    assert len(v[1:]) == len(seq), "LEN PROBLEM " +\
            str(len(v[1:])) + " " + str(len(seq))
   
    d=1
    for q, b in zip(v[1:],list(seq)):
        #only universal state has a self loop.
        if prev != q or q == laststate: 
            prob += A[prev, q]
            prob += D[prev, d]
            d = 1
        else:
            d += 1
        
        #if (q,b) not in B:
            #    print q, b, d, "not found"
            #else:
                #print q, b, d
        prob += B[q, b]
        prev = q
    

    return prob


def print_dict(T):
    for k in sorted(T):
        if isinstance(T[k], float):
            print "%s:%4.3f" % (k, T[k])
        else:
            print "%s:%s" % (k, T[k])
    print

def read_alphabet(afile):
    S = []
    f = open(afile, "rU")
    #assume space separated words
    for line in f:
        for word in line.strip().split():
            S.append(word)
    f.close()
    return S

class OptionParser (optparse.OptionParser):
    def check_required (self, *opt):
        for op in opt:
            option = self.get_option(op)
            # Assumes the option's 'default' is set to None!
            if getattr(self.values, option.dest) is None:
                self.error("%s option not supplied" % option)

parser = OptionParser()
# help message is automatically provided
parser.add_option('-A', '--AfromFile', dest='AfromFile',
                    help='read alphabet from given file') 
parser.add_option('-a', '--alphabet', dest='alphabet',
                    help='which alphabet to use (protein, kvogue, dna)') 
parser.add_option('-f', '--format', dest='format',
                    help='file format string (e.g., fasta)') 
parser.add_option('-o', '--outtrace', action='store_true', dest='outtrace',
                    help='print viterbi states')
parser.add_option('-q', '--qseqfile', dest='qseqfile',
                    help='compute prob of bouchra\'s qseq file')
parser.add_option('-r', '--trainfile', dest='trainfile',
                    help='train file')
parser.add_option('-R', '--hmmorder', type='int', dest='hmmorder', 
                    help='order of HMM (#chars/symbol)') 
parser.add_option('-t', '--testfile', dest='testfile',
                    help='test file')
parser.add_option('-v', '--voguefile', dest='voguefile',
                    help='VOGUE model file (.vgh)')

options, args = parser.parse_args(sys.argv[1:])
parser.check_required('-v', '-t')

if options.hmmorder is None:
    options.hmmorder = 1 #default one char/symbol

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
    Sigma = read_alphabet (options.AfromFile)
    options.format = "spaced"

psyco.bind(HMM.viterbi)

t0 = time.time()
hmm = HMM(options.voguefile, Sigma)
t1 = time.time()
#print "A:"
#print_dict(hmm.A)
#print "B:"
#print_dict(hmm.B)
#print "D:"
#print_dict(hmm.D)
#print "maxgap:", hmm.maxgap
#print "load time", t1-t0

#hmm.write_model()
if options.qseqfile is not None:
    f = open(options.qseqfile, "rU")
    while(1):
        inline = f.readline().strip().split()
        if inline[0] == 'LOG':
            slen = inline[3]
            qseq = [int(x) for x in inline[4:]]
            qseq.insert(0,0)
            print slen, len(qseq)
            break

testfile = open(options.testfile, "rU")

PS = ParseSeq.ParseSeq(testfile, options.format, options.hmmorder, kvogue)
PS.parse_next_sequence()
while PS.seq:
    t0 = time.time()
    #print PS.id
    if options.qseqfile is None:
        v = hmm.viterbi(PS.seq)
        qseq = v[0]
        hscore = v[1]
        t1 = time.time()
    
    if options.outtrace:
        print hscore, qseq
    
    nullscoredefault = HmmerNull.get_score_default(list(PS.seq), Sigma)
    nullscore = nullscoredefault
    #if options.hmmorder == 1 and not kvogue and options.alphabet != "dna":
    #    nullscore = HmmerNull.get_score(list(PS.seq))
    vscore = verify_score(qseq, hmm.A, hmm.B, hmm.D, PS.seq, hmm.n-1)
    if options.qseqfile is not None:
        hscore = vscore

    print PS.id, hscore, nullscore, nullscoredefault,\
            hscore-nullscore, hscore-nullscoredefault, vscore,\
            t1-t0, len(PS.seq)
    PS.parse_next_sequence()


