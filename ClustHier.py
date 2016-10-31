from BinaryTree import BinaryTree

class ClustHier:
    
    def __init__(self, matdists):
        self.matdists = matdists
    
    def distance (self, tree1, tree2):
        c1 = tree1.getCluster()
        c2 = tree2.getCluster()
        sd = 0.0
        for i in range(len(c1)):
            for j in range(len(c2)):
                sd += self.matdists.getValue(c1[i], c2[j])
        return sd/(len(c1)*len(c2)) 
    
    def executeClustering(self):
        trees = [ ]
        tableDist = self.matdists.copy()
        for i in range(self.matdists.numRows()):
            t = BinaryTree(i)
            trees.append(t)
        for k in range(self.matdists.numRows(), 1, -1):
            mins = tableDist.minDistIndexes()
            i,j = mins[0] , mins[1]
            n = BinaryTree(-1, tableDist.getValue(i, j)/2.0, trees[i], trees[j])
            if k>2:
                trees.pop(i)
                trees.pop(j)
                tableDist.removeRow(i)
                tableDist.removeRow(j)
                tableDist.removeCol(i)
                tableDist.removeCol(j)
                dists = [ ]
                for x in range(len(trees)): dists.append(self.distance(n, trees[x]))
                tableDist.addRow(dists)
                cdists = [ ]
                for y in range(len(dists)): cdists.append(dists[y])
                cdists.append(0.0)
                tableDist.addCol(cdists)
                trees.append(n)
            else: return n


