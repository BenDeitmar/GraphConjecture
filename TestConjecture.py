import numpy as np
from math import factorial
from itertools import product

def doubleFactorial(m):
    if m<-1:
        raise ValueError('doubleFactorial is not defined for {}'.format(m))
    elif m == -1:
        return 1
    elif m==0 or m==1:
        return 1
    else:
        return m*doubleFactorial(m-2)


def AdjacencyMatrixFromEdges(Edges,n):
    A = np.zeros((n,n))
    for e in Edges:
        A[e[0]-1,e[1]-1] += 1
        A[e[1]-1,e[0]-1] += 1
    for i in range(n):
        A[i][i] = A[i][i]/2
    return A

#recursively finds the set of vertices which are connected to the root
def recConnectedVertices(EdgeSet,root,initialSet):
    if EdgeSet:
        S = initialSet.copy()
        ToBeDiscarded=[]
        for e in EdgeSet:
            if e[0] in S:
                S.add(e[1])
                ToBeDiscarded.append(e)
            elif e[1] in S:
                S.add(e[0])
                ToBeDiscarded.append(e)
        if ToBeDiscarded:
            NewEdgeSet = EdgeSet.copy()
            for e in ToBeDiscarded:
                NewEdgeSet.discard(e)
            return recConnectedVertices(NewEdgeSet,root,S)
        else:
            return S
    else:
        return initialSet

def isConnected(EdgeSet,n):
    return len(recConnectedVertices(EdgeSet,1,set([1])))==n


class IgnorantUndirectedMultiGraph():
    def __init__(self, n, N, AdjacencyMatrix):
        self.n = n
        self.N = N
        self.A = AdjacencyMatrix

    #returns AdjacencyMatrix when the last vertex ist fused with vertex m
    def fused(self,m):
        n = self.n
        A = self.A.copy()
        for i in range(n):
            A[i][i] *= 2
        P = np.zeros((n-1,n))
        for i in range(n-1):
            P[i][i] = 1
        P[m][n-1] = 1

        A_new = np.matmul(P,A)
        A_new = np.matmul(A_new,P.transpose())
        for i in range(n-1):
            A_new[i][i] = A_new[i][i]//2
        return A_new

    def getLooplessDegree(self,i):
        return sum(self.A[i])-self.A[i,i]

    def getInclusiveDegree(self,i):
        return sum(self.A[i])+self.A[i,i]

    def getVertexDegree(self,i):
        return sum(self.A[i])

    def getLaplacian(self):
        n = self.n
        Q = np.multiply(self.A.copy(),-1)
        for i in range(n):
            Q[i][i] = self.getLooplessDegree(i)
        return Q

    def getSpanningTreeNumber(self):
        if self.n == 1:
            return 1
        Q = self.getLaplacian()
        #print(Q)
        Q = np.array(list(map(lambda L: L[:-1],Q[:-1])))
        return round(np.linalg.det(Q))

    def getWeight(self):
        n = self.n
        A = self.A
        W = 1
        for i in range(n):
            W *= factorial(self.getInclusiveDegree(i)-1)
            W *= doubleFactorial(2*int(A[i][i])-1)
        for j in range(1,n+1):
            W = W/factorial(2*A[j-1][j-1])
            for i in range(1,j):
                W = W/factorial(A[i-1][j-1])
        return W

    def __hash__(self):
        return hash((self.n,tuple(map(lambda L: tuple(L),self.A))))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.n == other.n and (self.A == other.A).all():
                return True
        return False


def checkFormulaOnGraph(G,P):
    n = G.n
    N = G.N
    leftHandSide = 2*(N+1-n)*G.getSpanningTreeNumber()*G.getWeight()
    rightHandSide = sum(map(lambda GP: GP[2]*GP[3],P))
    return (leftHandSide, rightHandSide)


def makeGeneology(n,N):
    ParentSet = set()
    L = list(map(lambda arg: product(range(1,n+2),range(1,n+2)), range(N)))
    for Edges in product(*L):
        A = AdjacencyMatrixFromEdges(Edges,n+1)
        G_Parent = IgnorantUndirectedMultiGraph(n+1,N,A)
        if not G_Parent in ParentSet and isConnected(set(Edges),n+1): #checking for connectednes is expensive
            #print(G_Parent.A)
            ParentSet.add(G_Parent)
    Geneology = dict()
    for G_Parent in ParentSet:
        t_Parent = G_Parent.getSpanningTreeNumber()
        W_Parent = G_Parent.getWeight()
        for m in range(n):
            G_Child = IgnorantUndirectedMultiGraph(n,N,G_Parent.fused(m))
            tup = (G_Parent,m,t_Parent,W_Parent)
            if G_Child in Geneology.keys():
                Geneology[G_Child].append(tup)
            else:
                Geneology[G_Child] = [tup]
    return Geneology

def CheckFormula(n,N):
    if n>N:
        raise Exception("Error: There are no connected parents for n>N")
    Geneology = makeGeneology(n,N)
    AnyExceptions = False
    for G_Ghild in Geneology.keys():
        print("Checking graph:")
        print(G_Ghild.A)
        l,r = checkFormulaOnGraph(G_Ghild,Geneology[G_Ghild])
        print(l==r, l,'=',r)
        if l!=r:
            AnyExceptions = True
    if AnyExceptions:
        print("Exception found!!!")
    else:
        print("No Exception")

def FindConnection(n,N):
    if n>N:
        raise Exception("Error: There are no connected parents for n>N")
    Geneology = makeGeneology(n,N)
    for G_Ghild in Geneology.keys():
        l,r = checkFormulaOnGraph(G_Ghild,Geneology[G_Ghild])
        print(G_Ghild.A,l/r)



if __name__ == "__main__":

    max_N = 4
    for N in range(1,max_N+1):
        for n in range(1,N+1):
            print("########################################")
            print("n={} N={}".format(n,N))
            CheckFormula(n,N)
