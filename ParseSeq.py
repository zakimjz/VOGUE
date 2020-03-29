import HmmerNull
from Bio import SeqIO


class ParseSeq:
    def __init__(self, handle, format, hmmorder, kvogue):
        self.id = None
        self.seq = []
        self._handle = handle
        self._format = format
        self._hmmorder = hmmorder
        self._kvogue = kvogue
        if self._format == "fasta":
            self._record = SeqIO.parse(self._handle, "fasta")
        elif self._format == "spaced":
            # read all lines into local mem
            self._lines = []
            self._ids = []
            for line in self._handle:
                if line[0] != '>':
                    A = [w for w in line.strip().split()]
                    self._lines.append(A)
                else:
                    id = line.strip().split()[1]
                    self._ids.append(id)
            self._curidx = 0  # init to first line

    def parse_next_sequence(self):
        self.id = None
        self.seq = []
        if self._format == "fasta":
            try:
                record = next(self._record)
            except StopIteration:
                return
            #print record.id
            self.id = record.id
            for i in range(0, len(record.seq)-self._hmmorder+1,
                           self._hmmorder):
                # wd is a word of len hmmorder, counts as one sym
                wd = ''
                for w in record.seq[i:i+self._hmmorder]:
                    if self._kvogue:
                        w = HmmerNull.KREP[w]
                    wd += str(w)
                self.seq.append(wd)
        elif self._format == "spaced":
            if self._curidx >= len(self._lines):
                return
            else:
                self.id = self._ids[self._curidx]
                self.seq = self._lines[self._curidx]
                self._curidx += 1
        #print self.id, self.seq
