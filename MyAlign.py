
class MyAlign:
    def __init__(self, lseqs, tipo="protein"):
        self.listseqs = lseqs
        self.tipo = tipo

    def __len__(self):
        return len(self.listseqs[0])
 
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2:
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int: 
            return self.listseqs[n]
 
    def __str__(self):
        res = ""
        for seq in self.listseqs: res += "\n" + seq
        return res 

    def consensus (self):
        cons = ""
        for i in range(len(self)): # cada coluna
            cont = { }
            for k in range(len(self.listseqs)): # linha / seq
                c = self.listseqs[k][i]
                if c in cont: 
                    cont[c] = cont[c] + 1
                else: 
                    cont[c] = 1
            maximum = 0
            cmax = None
            for ke in cont.keys():
                if ke != "-" and cont[ke] > maximum:
                    maximum = cont[ke]
                    cmax = ke
            cons = cons + cmax
        return cons 
        
    def column(self,i): 
        res = [ ]
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][i])
        return res


if __name__ == "__main__":
    alig = MyAlign(["ATGA-A","AA-AT-"], "dna")
    print (alig)
    print (len(alig))
    print (alig[1,1])
    print (alig[0,2])
    print (alig.consensus())
    
    