from MatrixNum import MatrixNum
from MySeq import MySeq
from SubstMatrix import SubstMatrix
from AlignSeq import AlignSeq
from ClustHier import ClustHier

class UPGMA:
    
    def __init__(self, seqs, alseq):
        self.seqs = seqs
        self.alseq = alseq
        self.matdist = MatrixNum(len(seqs), len(seqs))
        self.criaMatDists() 
        
    def criaMatDists(self):
        for i in range(len(self.seqs)):
            self.matdist.setValue(i, i, 0.0)
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.needlemanWunsch(s1, s2)
                alin = self.alseq.recoverAlignment()
                ncd = 0
                for k in range(len(alin)):
                    col = alin.column(k)
                    if (col[0] != col[1]): ncd += 1
                self.matdist.setValue(i, j, ncd)
                self.matdist.setValue(j, i, ncd)         
        
    def run(self):
        ch = ClustHier(self.matdist)
        arv= ch.executeClustering()
        return arv
        
        
   
def teste():
    seq1 = MySeq("ATAGCGAT")
    seq2 = MySeq("ATAGGCCT")
    seq3 = MySeq("CTAGGCCC")
    sm = SubstMatrix()
    sm.createFromMatchPars(3, -1, "ACGT")
    alseq = AlignSeq(sm, -8)
    up = UPGMA([seq1, seq2, seq3], alseq)
    arv = up.run()
    arv.printtree() 


if __name__ == '__main__':
    teste()