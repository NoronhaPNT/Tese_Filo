

import math
import numpy

class MatrixNum:

    def __init__(self, rows, cols):
        self.mat = []
        for i in range(rows):
            self.mat.append([])
            for j in range(cols):
                self.mat[i].append(0.0)

    def __getitem__(self, n):
        return self.mat[n]
        
    def numRows (self):
        return len(self.mat)
    
    def numCols (self):
        return len(self.mat[0])
    
    def getValue (self, i, j):
        return self.mat[i][j]
    
    def setValue(self, i, j, value):
        self.mat[i][j] = value
    
    def printmat(self):
        for r in self.mat: print (r)

    def sum (self):
        s = 0.0;
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                s +=  self.mat[i][j]        
        return s

    def mean (self):
        return self.sum() / (self.numRows() * self.numCols())

    def maximum (self):
        m = self.mat[0][0]
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if (self.mat[i][j] > m):
                    m = self.mat[i][j]        
        return m
    
    def minimum (self):
        m = self.mat[0][0]
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if (self.mat[i][j] < m):
                    m = self.mat[i][j]        
        return m
   
    def maxIndexes (self):
        m = self.mat[0][0]
        res = (0,0)
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if (self.mat[i][j] > m):
                    m = self.mat[i][j]
                    res = (i,j)        
        return res
    
    def minIndexes (self):
        m = self.mat[0][0]
        res = (0,0)
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if (self.mat[i][j] < m):
                    m = self.mat[i][j]
                    res = (i,j)        
        return res
    
    def minDistIndexes (self):
        m = self.mat[1][0]
        res= (1,0)
        for i in range(self.numCols()):
            for j in range(i+1, self.numRows()):
                if self.mat[j][i] < m:
                    m = self.mat[j][i]
                    res = (j, i)
        return res
    
    def meanRow (self, li):
        s = 0
        for k in range(len(self.mat[li])):
            s += self.mat[li][k]
        return s / len(self.mat[li])

    def maxRow (self, li):
        m = self.mat[li][0]
        for k in range(1,len(self.mat[li])):
            if self.mat[li][k] > m:
                m = self.mat[li][k]
        return m

    def minRow (self, li):
        m = self.mat[li][0]
        for k in range(1,len(self.mat[li])):
            if self.mat[li][k] < m:
                m = self.mat[li][k]
        return m

    def meanRows(self):
        res = []
        for r in range(self.numRows()):
            res.append(self.meanRow(r))
        return res

    def minRows(self):
        res = []
        for r in range(self.numRows()):
            res.append(self.minRow(r))
        return res

    def maxRows(self):
        res = []
        for r in range(self.numRows()):
            res.append(self.maxRow(r))
        return res


    def meanCol(self, lc):
        s = 0
        for k in range(self.numCols()):
            s += self.mat[k][lc]
        return s / self.numCols()
    
    def maxCol(self, lc):
        m = self.mat[0][lc]
        for c in range(1, self.numCols()):
            if self.mat[c][lc] > m:
                m = self.mat[c][lc]
        return m
    
    def minCol(self, lc):
        m = self.mat[0][lc]
        for c in range(1, self.numCols()):
            if self.mat[c][lc] < m:
                m = self.mat[c][lc]
        return m
    
    def meanCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.meanCol(c))
        return res

    def maxCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.maxCol(c))
        return res
    
    def minCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.minCol(c))
        return res
    
    def applyLog(self, b):
        newm = MatrixNum(self.numRows(), self.numCols())
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                newm[i][j] = math.log(self.getValue(i, j)) / math.log(b)
        return newm       

    def normalize(self):
        sd = numpy.std(self.mat)
        m = self.mean()
        newm = MatrixNum(self.numRows(), self.numCols())
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                newm[i][j] = (self.getValue(i, j) - m) / sd
        return newm    
       
    def addRow(self, newrow):
        self.mat.append(newrow)

    def addCol(self, newcol):
        for r in range(self.numRows()):
            self.mat[r].append(newcol[r])

    def removeRow(self, ind):
        del self.mat[ind]

    def removeCol(self, ind):
        for r in range(self.numRows()):
            del self.mat[r][ind]

    def copy(self):
        newm = MatrixNum(self.numRows(), self.numCols())
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                newm.mat[i][j] = self.mat[i][j]
        return newm
    
def test():
    m = MatrixNum(3,2)
    m.setValue(2,1,3.0)
    m.setValue(1,0,2.0)
    m.setValue(0,1,5.0)
    m.setValue(1,1,-2.0)
    m.printmat()
    print (str(m.mean()))
    print (str(m.maximum()))
    print (str(m.minIndexes()))
    print (str(m.meanRows()))
    print (str(m.meanCols()))
    m.removeRow(1)
    m.printmat()
    m.addRow([2.0,4.0])
    m.printmat()
    nm = m.copy()
    nm.printmat()
    nm.addCol([3.0,4.0,5.0])
    nm.printmat()
    nm.removeCol(0)
    nm.printmat()
    print (str(nm.meanRows()))

if __name__ == '__main__': # testa se esta script e lancada
    test()