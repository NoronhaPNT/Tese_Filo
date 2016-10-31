
import sys

class BinaryTree:

    def __init__(self, val, dist=0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right

        
    def getCluster(self):
        res = []
        if self.value >= 0:
            res.append(self.value)
        else:
            if (self.left != None):
                res.extend(self.left.getCluster())
            if (self.right != None): 
                res.extend(self.right.getCluster())
        return res
    
    def printtree(self):
        self.printtreerec(0, "Root")
    
    def printtreerec (self, level, side):
        for i in range(level): sys.stdout.write("\t")
        al = self.getCluster();
        sys.stdout.write(side + ":" + str(al)+ " Dist.: " + str(self.distance) + "\n")
        if self.value < 0:
            if (self.left != None): 
                self.left.printtreerec(level+1, "Left")
            else: 
                sys.stdout.write("Null")
            if (self.right != None): 
                self.right.printtreerec(level+1, "Right")
            else: 
                sys.stdout.write("Null\n")

    def size(self):
        numleafes = 0
        numnodes = 0
        if self.value >= 0:
            numleafes = 1
        else: 
            if (self.left != None):
                resl = self.left.size()
            else: resl = (0,0)
            if (self.right != None):  
                resr = self.right.size() 
            else: resr = (0,0)
            numnodes += (resl[0] + resr[0] + 1)
            numleafes += (resl[1] + resr[1])
        return numnodes, numleafes

    def existsLeaf(self, leafnum):
        if self.value >= 0:
            if self.value == leafnum:
                return True
            else: return False
        else:
           if self.left != None:
               resl = self.left.existsLeaf(leafnum)
               if resl == True: return True
           if self.right != None:
               resr = self.right.existsLeaf(leafnum) 
               if resr == True: return True
        return False
        
    def commonAncestor(self, leaf1, leaf2):
        if self.value >= 0: return None
        if self.left.existsLeaf(leaf1):
            if self.left.existsLeaf(leaf2):
                return self.left.commonAncestor(leaf1, leaf2)
            if self.right.existsLeaf(leaf2):
                return self
            return None       
        if self.right.existsLeaf(leaf1):
            if self.right.existsLeaf(leaf2):
                return self.right.commonAncestor(leaf1, leaf2)
            if self.left.existsLeaf(leaf2):
                return self
        return None

def test():              
    a = BinaryTree(1)
    b = BinaryTree(2)
    c = BinaryTree(3)
    d = BinaryTree(4)
    e = BinaryTree(-1, 2.0, a, b)
    f = BinaryTree(-1, 3.0, c, d)
    g = BinaryTree(-1, 4.0, e, f)
    g.printtree()
    print (g.size())
    print (g.existsLeaf(3))
    print (g.existsLeaf(7))
    g.commonAncestor(1,2).printtree()
    g.commonAncestor(1,4).printtree()
    print (g.commonAncestor(5,1))

if __name__ == '__main__':
    test()
    
# Result for test
# Root:[1, 2, 3, 4] Dist.: 4.0
#     Left:[1, 2] Dist.: 2.0
#         Left:[1] Dist.: 0
#         Right:[2] Dist.: 0
#     Right:[3, 4] Dist.: 3.0
#         Left:[3] Dist.: 0
#         Right:[4] Dist.: 0