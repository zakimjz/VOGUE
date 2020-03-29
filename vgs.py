import math
import ParseSeq

def compute_interleaved_pos(posinfo, key1, key2):
    diff = 0.0
    #build a map for key1, px
    pidict = {}
    for (x,y) in posinfo[key1]:
        if x not in pidict: pidict[x] = {}
        pidict[x][y] = True

    #check for i,j in key2 poslist if x,i and j,y exist in key1
    for (i,j) in posinfo[key2]:
        if j in pidict: #if key j exists than j,y exists
            found = False
            for x in pidict: #check if x,i exists
                if i in pidict[x]:
                    found = True
                    break
            if not found: #pair not found, increase diff
                diff += 1
        else: diff += 1

    return diff

def remove_artifact_seqs(posinfo, thresh):
    if thresh is None:
        thresh = 0.2 #consider artifact if counts differ by <= thresh
    keys = posinfo.keys();
    for key1 in keys:
        key2 = (key1[1], key1[0])
        if key1 == key2 or key1 not in posinfo or key2 not in posinfo:
            continue
        sup1 = len(posinfo[key1])
        sup2 = len(posinfo[key2])
        if sup1 >= sup2: 
            rkey = key2 #remove the smaller support key
            diff = compute_interleaved_pos(posinfo, key1, key2)
        else: 
            rkey = key1 #remove the smaller support key
            diff = compute_interleaved_pos(posinfo, key2, key1)
        
        ratio = diff/len(posinfo[rkey]) 

        if ratio <= thresh:
            #print "P1", posinfo[key1]
            #print "P2", posinfo[key2]
            del posinfo[rkey]
            print "DEL KEY", rkey, diff, ratio


def vgs(infile, maxgap, minsup, minsupratio, wminsup,
        format, kvogue, hmmorder, artifacts, trainsize=None):
    itemcount = {} #hold 1-seq freq
    seqcount = {} #hold 2-seq freq
   
    #find frequent 2-seq
    seqcnt = 0
    handle = open(infile, "rU")
    PS = ParseSeq.ParseSeq(handle, format, hmmorder, kvogue)
    PS.parse_next_sequence()
    while PS.seq:
        seqcnt += 1
        if trainsize is not None and seqcnt > trainsize: break
        for x in PS.seq:
            if x not in itemcount:
                itemcount[x] = 0.0
            itemcount[x] += 1

        for i in range(len(PS.seq)-1):
            for j in range(i+1, min(len(PS.seq), i+maxgap+2)):
                pair = (PS.seq[i], PS.seq[j])
                #print pair
                if not seqcount.has_key(pair):
                    seqcount[pair] = [None]

                if seqcount[pair][-1] != PS.id:
                    seqcount[pair].append(PS.id)
        PS.parse_next_sequence()
                    
    handle.close()

    #print seqcount
    #delete infreq pairs
    if minsupratio is not None:
        minsup = math.ceil(minsupratio*seqcnt)

    print "MINSUP: ", minsup, seqcnt
    keys = seqcount.keys();
    for key in keys:
        sup = len(seqcount[key])-1
        if sup < minsup:
            del seqcount[key]

    posinfo = {} #info about the positions (for artifact removal)
    gapinfo = {} #info about the gaps
    handle = open(infile, "rU")
    seqcnt = 0
    PS = ParseSeq.ParseSeq(handle, format, hmmorder, kvogue)
    PS.parse_next_sequence()
    while PS.seq:
        seqcnt += 1
        if trainsize is not None and seqcnt > trainsize: break

        for i in range(len(PS.seq)-1):
            for j in range(i+1, min(len(PS.seq), i+maxgap+2)):
                frag = PS.seq[i:j+1]
                #print i, j, frag
                gap = j-i-1;
                key = frag[0],frag[-1]

                if key in seqcount: #known to be frequent
                    if artifacts:
                        if key not in posinfo: posinfo[key] = []
                        posinfo[key].append((i,j))

                    if key not in gapinfo: #initialize
                        #[0] has count, [1] count by gap, [2] gap symbols
                        gapinfo[key] = [0.0, [0.0]*(maxgap+1), {}] 

                    gapinfo[key][0] += 1 #increment weighted sup
                    gapinfo[key][1][gap] += 1
                    
                    for g in range(1,gap+1):
                        gapsym = PS.seq[i+g]
                        if gapsym not in gapinfo[key][2]:
                            gapinfo[key][2][gapsym] = 0.0
                        gapinfo[key][2][gapsym] += 1
        PS.parse_next_sequence()
                        
    handle.close()
    
    if artifacts: #C-VOGUE
        remove_artifact_seqs(posinfo, artifacts)

    keys = seqcount.keys()
    for key in keys: #known to be freq
        if artifacts and key not in posinfo:
            del gapinfo[key]

        if wminsup is not None: #apply threshold on weighted minsup
            if gapinfo[key][0] < wminsup:
                del gapinfo[key]


    return itemcount, gapinfo
