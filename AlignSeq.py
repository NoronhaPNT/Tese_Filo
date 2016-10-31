
from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix

class AlignSeq:

    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.S = None
        self.T = None
        self.seq1 = None
        self.seq2 = None
        
    def scorePos (self, c1, c2):
        if c1 == "-" or c2=="-":
            return self.g
        else:
            return self.sm[c1,c2]
        
    def scoreAlin (self, alin):
        res = 0;
        for i in range(len(alin)):
            res += self.scorePos (alin[0][i], alin[1][i])
        return res
    
    def needlemanWunsch (self, seq1, seq2):
        if (seq1.tipo != seq2.tipo): return None
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        for j in range(1, len(seq2)+1):
            self.S[0].append(self.g * j)
            self.T[0].append(3)
        for i in range(1, len(seq1)+1):
            self.S.append([self.g * i])
            self.T.append([2])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i][j] + self.scorePos (seq1[i], seq2[j])
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                self.S[i+1].append(max(s1, s2, s3))
                self.T[i+1].append(max3t(s1, s2, s3))
        return self.S[len(seq1)][len(seq2)]
    
    def recoverAlignment (self):
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)
        while i>0 or j>0:
            if self.T[i][j]==1:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = self.seq2[j-1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.seq2[j-1] + res[1] 
                j -= 1
            else:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
        return MyAlign(res, self.seq1.tipo)
    
    def smithWaterman (self, seq1, seq2):
        if (seq1.tipo != seq2.tipo): return None
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        maxscore = 0
        for j in range(1, len(seq2)+1):
            self.S[0].append(0)
            self.T[0].append(0)
        for i in range(1, len(seq1)+1):
            self.S.append([0])
            self.T.append([0])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i][j] + self.scorePos(seq1[i], seq2[j]) 
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                b = max(s1, s2, s3)
                if b <= 0:
                    self.S[i+1].append(0)
                    self.T[i+1].append(0)
                else:
                    self.S[i+1].append(b)
                    self.T[i+1].append(max3t(s1, s2, s3))
                    if b > maxscore: 
                        maxscore = b
        return maxscore

    def recoverLocalAlignment (self):
        res = ["", ""]
        maxscore = 0
        maxrow = 0
        maxcol = 0
        for i in range(1,len(self.S)):
            for j in range(1, len(self.S[i])):
                if self.S[i][j] > maxscore:
                    maxscore = self.S[i][j]
                    maxrow = i
                    maxcol = j
        i = maxrow
        j = maxcol
        while i>0 or j>0:
            if self.T[i][j]==1:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = self.seq2[j-1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                res[0] = "-" + res[0];
                res[1] = self.seq2[j-1] + res[1]; 
                j -= 1
            elif self.T[i][j] == 2:
                res[0] = self.seq1[i-1] + res[0];
                res[1] = "-" + res[1]; 
                i -= 1
            else: break
        return MyAlign(res, self.seq1.tipo)

    
    def alinhaSeqComConjunto(self, query, conjunto):
        scores = []
        res = [[0,0]]*len(conjunto)
        k = 0
        for s in conjunto:
            scores[k] = self.smithWaterman(query, s)
            k = k + 1
        for i in range(len(conjunto)):
            j = 0;
            while j < i and scores[ res[j][0] ] > scores[i]: 
                j = j + 1
            for k in range(i,j,-1):
                res[k][0] = res[k-1][0]
            res[j][0] = i
        for i in range(len(conjunto)):
            res[i][1] = scores[ res[i][0] ]
        return res    
    
    def alinhaConjuntos(self, c1, c2):
        res = []
        for i in range(len(c1)):
            maxscore = -1
            bestseq = None
            for j in range(len(c2)):
                score = self.smithWaterman(c1[i], c2[j])
                if score > maxscore:
                    maxscore = score
                    bestseq = c2[j]
            res.append(bestseq)
        return res
    
    def alinhaProtDNA(self, prot, dna):
        res = None
        orfs = dna.orfs();
        prots = []
        for i in range(len(orfs)): 
            prots.append( orfs[i].extraiTodasProts() )
        mx = 0
        for p in prots:
            score = self.smithWaterman(prot.seq, p.seq)
            if score > mx:
                mx = score
                res = self.getLocalAlignment()
        return res;
   
    
def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3
    
def printMat (mat):
    for i in range(0, len(mat)):
        print (mat[i])
        

"""
def test1():
    seq1 = MySeq("NCRD","protein")
    seq2 = MySeq("CRNC","protein")
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    alin = AlignSeq(sm, -4)
    print alin.needlemanWunsch(seq1, seq2)
    printMat(alin.S)
    print
    printMat(alin.T)
    print alin.recoverAlignment()
"""
# Result from test 1
# 7
# [0, -4, -8, -12, -16]
# [-4, -3, -4, -2, -6]
# [-8, 5, 1, -3, 7]
# [-12, 1, 10, 6, 3]
# [-16, -3, 6, 11, 7]
# 
# [0, 3, 3, 3, 3]
# [2, 1, 1, 1, 3]
# [2, 1, 3, 3, 1]
# [2, 2, 1, 3, 2]
# [2, 2, 2, 1, 3]
# 
# NCRD-
# -CRNC
"""
def test2():
    seq1 = MySeq("ATGATATGATGATT")
    seq2 = MySeq("GATGAATAGATGTGT")
    sm = SubstMatrix()
    sm.createFromMatchPars(3, -1, "ACGT")
    alin = AlignSeq(sm, -3)
    print alin.smithWaterman(seq1, seq2)
    printMat(alin.S)
    print alin.recoverLocalAlignment()
"""
# Results from test 2
# 25
# [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# [0, 0, 3, 0, 0, 3, 3, 0, 3, 0, 3, 0, 0, 0, 0, 0]
# [0, 0, 0, 6, 3, 0, 2, 6, 3, 2, 0, 6, 3, 3, 0, 3]
# [0, 3, 0, 3, 9, 6, 3, 3, 5, 6, 3, 3, 9, 6, 6, 3]
# [0, 0, 6, 3, 6, 12, 9, 6, 6, 4, 9, 6, 6, 8, 5, 5]
# [0, 0, 3, 9, 6, 9, 11, 12, 9, 6, 6, 12, 9, 9, 7, 8]
# [0, 0, 3, 6, 8, 9, 12, 10, 15, 12, 9, 9, 11, 8, 8, 6]
# [0, 0, 0, 6, 5, 7, 9, 15, 12, 14, 11, 12, 9, 14, 11, 11]
# [0, 3, 0, 3, 9, 6, 6, 12, 14, 15, 13, 10, 15, 12, 17, 14]
# [0, 0, 6, 3, 6, 12, 9, 9, 15, 13, 18, 15, 12, 14, 14, 16]
# [0, 0, 3, 9, 6, 9, 11, 12, 12, 14, 15, 21, 18, 15, 13, 17]
# [0, 3, 0, 6, 12, 9, 8, 10, 11, 15, 13, 18, 24, 21, 18, 15]
# [0, 0, 6, 3, 9, 15, 12, 9, 13, 12, 18, 15, 21, 23, 20, 17]
# [0, 0, 3, 9, 6, 12, 14, 15, 12, 12, 15, 21, 18, 24, 22, 23]
# [0, 0, 0, 6, 8, 9, 11, 17, 14, 11, 12, 18, 20, 21, 23, 25]
# 
# ATGATAT-GATGATT
# ATGA-ATAGATGTGT
    

def lcs (seq1, seq2):
    sm = SubstMatrix()
    sm.createFromMatchPars(1,0, seq1.alphabet())
    aseq = AlignSeq(sm, 0)
    aseq.needlemanWunsch(seq1, seq2)
    alin = aseq.recoverAlignment ()
    sizeal = len(alin[0])
    lcs = ""
    for i in range(sizeal):
        if alin[0][i] != "-" and alin[0][i] == alin[1][i]:
            lcs += alin[0][i]
    return lcs

def testLCS():
    s1 = MySeq("ATTAGCT")
    s2 = MySeq("ATAAGCT")
    print (lcs(s1,s2))
    
def edit_distance(seq1, seq2):
    sm = SubstMatrix()
    sm.createFromMatchPars(0,-1,seq1.alphabet())
    aseq = AlignSeq(sm, -1)
    sc = aseq.needlemanWunsch(seq1, seq2)
    return -sc

def testED():
    s1 = MySeq("ATTAGCT")
    s2 = MySeq("ATTAAAGCT")
    print (edit_distance(s1,s2))

if __name__ == "__main__":   
    test1()
    test2()