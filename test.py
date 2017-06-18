import random
import math
import networkx as nx
import igraph as ig
import itertools
import numpy as np
import timeit
import bisect
from BinPackerDynamic import Binpacker,Item
from kernighanLin import kernighan_lin_bisection
from Tree import Tree,TreeUtil

import matplotlib.pyplot as plt

class GraphPartitioning(object):

    def __init__(self,epsilon,n,k,G):
        self.epsilon = epsilon
        self.n = n
        self.k = k
        self.G = G
        self.feasibleSets = set()
        self.unfeasibleSets = set()

        print("Is it connected?:   " + str(nx.is_connected(G)))
        if (nx.is_connected(G) == False):
            print("Removing Singols....")
            listOfNodes = G.nodes()
            for x in range(len(listOfNodes)):
                if (G.degree(str(listOfNodes[x])) == 0):
                    G.remove_node(str(listOfNodes[x]))
            print("Is it connected now?:   " + str(nx.is_connected(G)))

    def fromNetworkXtoiGraph(self, subTree):  # Conversion from networkX data structure to igraph data
        newlist = []
        weightVector = []
        edgeList = subTree.data.edges()
        for x in range(len(edgeList)):  # Conversion from NetworkX Graph to igraph data structure
            couple = []
            couple.append(
                int(edgeList[x][0]))  # the edges has strucure [(n1,n2,{}),....] {} indicates addition attribute
            couple.append(int(edgeList[x][1]))
            newlist.append(couple)
        return ig.Graph(len(subTree.data.nodes()), newlist)

    def createSubTree(self, partition, subTree, threshold, G):  #
        X = G.subgraph(str(x) for x in partition[0])  # spotting the two sub-sets divided inside the original Graph G
        Y = G.subgraph(str(x) for x in partition[1])
        leftChild = Tree()  # declaration offsprings
        leftChild.data = X  # Right Child data assignment
        leftChild.left = None
        leftChild.right = None
        rightChild = Tree()
        rightChild.data = Y  # Left Child data assignment
        rightChild.left = None
        rightChild.right = None
        subTree.left = leftChild  # childreen relation assignment
        subTree.right = rightChild
        if (len(partition[0]) > threshold):  # progressive recursive approach from root to leaves
            self.division(leftChild, threshold, G)
        if (len(partition[
                    1]) > threshold):  # treshold is evaluated in order to avoid the creation of parts that need to be removed  (epsilon*n)/ 3*k
            self.division(rightChild, threshold, G)

    def division(self, subTree, threshold,G):
        if (subTree.data.number_of_nodes() <= 1):  # if subtree is a leaf, return it
            return subTree
        else:
            balanceCut = int(round((((self.epsilon / 3) / (
            1 + self.epsilon / 3)) / 2) * subTree.data.number_of_nodes()))  # Evaluation of epsilon' according to formula provided in paper (Section "3. Algortithm and Analysis")
            listOfNodes = subTree.data.nodes()  # Selection of two randon nodes in order to divide the graph in two separate parts
            nodeSource = random.choice(listOfNodes)  # source node
            listOfNodes.remove(nodeSource)
            nodeDestination = random.choice(listOfNodes)  # destination node
            graphNew = self.fromNetworkXtoiGraph(subTree)
            mc = graphNew.mincut(int(nodeSource), int(nodeDestination))  # cutting method of igraph library
            partition = mc.partition
            while True:  # igraph has a continus declaration of nodes, so it is not possible to declare a non-continue set of numbers
                cont = 0  # ex: set(3,4) -> creates the set (1,2,3,4) [1,2] needs to be removed
                n = 0
                for x in range(len(partition[0])):  # This "for" fix the problem of igraph continuity (2 line before)
                    if (str(partition[0][x - n]) in subTree.data.nodes()):
                        cont += 1
                    else:
                        del partition[0][x - n]  # removing the additional elements
                        n += 1
                if (cont >= balanceCut and len(
                        partition[1]) >= balanceCut):  # if the two subsets respect the bipartion constraints of
                    break  # (e')/2*|V| we have created the two subsets.
                else:
                    listOfNodes = subTree.data.nodes()  # otherwise we search another bipartion, according the
                    nodeSource = random.choice(listOfNodes)  # constraints
                    listOfNodes.remove(nodeSource)
                    nodeDestination = random.choice(listOfNodes)
                    mc = graphNew.mincut(int(nodeSource), int(nodeDestination))  # defining a new partition
                    partition = mc.partition
            self.createSubTree(partition, subTree, threshold, G)

    def bilancedPartion(self,subTree, threshold, G):
        monoWeight = True
        if (subTree.data.number_of_nodes() <= 1):  # if subtree is a leaf, return it
            return subTree
        else:
            balanceCut = int(round(((self.epsilon / 3) / (
            1 + self.epsilon / 3)) / 2 * subTree.data.number_of_nodes()))  # Evaluation of epsilon according to formula provided in paper
            listOfNodes = subTree.data.nodes()  # Selection of two randon nodes in order to divide the graph in two separate parts
            print("Avvio del taglio")
            partition = kernighan_lin_bisection(subTree.data)
            X = G.subgraph(partition[0])  # spotting the two sub-sets of division inside the original Graph G
            Y = G.subgraph(partition[1])
            print(str(len(subTree.data)) + "  " + str(len(X)) + "  " + str(len(Y)))
            leftChild = Tree()  # declaration of offsprings
            leftChild.data = X  # Right Child data assignment
            leftChild.left = None
            leftChild.right = None
            rightChild = Tree()
            rightChild.data = Y  # Left Child data
            rightChild.left = None
            rightChild.right = None
            subTree.left = leftChild
            subTree.right = rightChild
            if (len(partition[0]) > threshold):
                self.bilancedPartion(leftChild, threshold,G)
            if (len(partition[1]) > threshold):
                self.bilancedPartion(rightChild, threshold,G)


    def removeBigNodes(self,root,threshold,listTrees):
        if(root.data.number_of_nodes()<threshold):
            listTrees.append(root)
            return
        else:
            self.removeBigNodes(root.left,threshold,listTrees)
            self.removeBigNodes(root.right, threshold, listTrees)

    def pruneTree(self,root):
        listTrees = []
        self.removeBigNodes(root,(1+self.epsilon)*(self.n/self.k),listTrees)
        return listTrees

    def getGVector(self):
        vector_gValue = []
        val = math.pow((1 + self.epsilon/2),
                       math.floor(math.log(self.epsilon*self.n/(3*self.k),
                                           (1+self.epsilon/2))))
        while(val < ((1+self.epsilon)*(self.n/self.k))):
            vector_gValue.append(val)
            val *= (1 + self.epsilon / 2)
        return vector_gValue

    def setFixedPartitions(self,listTrees,packer,GVector):
        vector_g = [0] * len(GVector)
        for tree in listTrees:
            if(tree.left == None and tree.right == None):
                numberOfNodes = tree.data.number_of_nodes()
                if(numberOfNodes<=GVector[0]):
                    packer.items.append(Item('A', round(GVector[0])))
                for i in range(len(GVector)-1):
                    if (numberOfNodes < GVector[i+1] and numberOfNodes >= GVector[i]):
                        packer.items.append(Item('A', round(GVector[i+1])))

    def calculateCost(self,graphData,left,right):            #calculate cost of divide two sets
        sum = 0
        leftNodes = left.data.node
        rightNodes = right.data.node
        smallPart = None
        bigPartPart = None
        if(len(leftNodes) > len(rightNodes)):
            smallPart = rightNodes
            bigPartPart = leftNodes
        else:
            smallPart = leftNodes
            bigPartPart = rightNodes
        for edge in graphData.edge:
            if(edge in smallPart):
                edges = graphData.edge[edge]
                for singEdge in edges:
                    if(singEdge in bigPartPart):
                        sum += edges[singEdge]['weight']
        return sum

    def calculateCostFinal(self,graphData,left,right):          #calculate final cost of a partitioning
        sum = 0
        leftNodes = left
        rightNodes = right
        smallPart = None
        bigPartPart = None
        if(len(leftNodes) > len(rightNodes)):
            smallPart = rightNodes
            bigPartPart = leftNodes
        else:
            smallPart = leftNodes
            bigPartPart = rightNodes
        for edge in graphData.edge:
            if(edge in smallPart):
                edges = graphData.edge[edge]
                for singEdge in edges:
                    if(singEdge in bigPartPart):
                        sum += edges[singEdge]['weight'] #= 1
        return sum

    def getPartitionInTree(self,listNodes,num,i):
        tempVett = []
        if(i[0]==num):
            temp = []
            for el in listNodes:
                temp.append(el.data.node.keys())
            return temp
        for k in range(len(listNodes)):
            if (listNodes[k].left != None):
                father = listNodes[k]
                listNodes.pop(k)
                listNodes.insert(k, father.left)
                listNodes.insert(k + 1, father.right)
                i[0] += 1
                tempVett = self.getPartitionInTree(listNodes, num, i)
                if(tempVett != None):
                    return tempVett
                listNodes.pop(k)
                listNodes.pop(k)
                listNodes.insert(k, father)

    def createPossiblePartitions(self,listNodes,listTreesCosts,cost):
        listTreesCosts.append(cost)
        for i in range(len(listNodes)):
            if(listNodes[i].left !=None):
                father = listNodes[i]
                listNodes.pop(i)
                listNodes.insert(i,father.left)
                listNodes.insert(i+1, father.right)
                newCost = self.calculateCost(father.data,father.left,father.right)
                self.createPossiblePartitions(listNodes,listTreesCosts,cost+newCost)
                listNodes.pop(i)
                listNodes.pop(i)
                listNodes.insert(i,father)

    def scanTrees(self,listTrees):
        listOfTotalCosts = []
        for tree in listTrees:
            temp = []
            listTreesCosts = []
            temp.append(tree)
            self.createPossiblePartitions(temp,listTreesCosts,0)
            listOfTotalCosts.append(listTreesCosts)
        return listOfTotalCosts

    def getPartitionCost(self,n,PossiblePartitionsCost):
        cost = 0
        for i in range(len(PossiblePartitionsCost)):
            module = len(PossiblePartitionsCost[-i - 1])
            cost += PossiblePartitionsCost[-i - 1][int(n % module)]
            n -= (n % module)
            n /= module
        return cost

    def createPartitionNodes(self,listTrees,n,PossiblePartitionsCost,GVector,packer,lastPartitions,indexLastPartitions):
        partition = []
        currentG = np.zeros(len(GVector))
        for i in range(len(listTrees)):#tree in listTrees[::-1]:
            temp = []
            temp.append(listTrees[-i-1])
            module = len(PossiblePartitionsCost[-i-1])
            tempPartition = self.getPartitionInTree(temp,n%module,[0])
            indexLastPartitions[i] = n
            lastPartitions [i] = tempPartition
            n -= (n%module)
            n /= module
            if(tempPartition != None):
                for el in tempPartition:
                    partition.append(el)
                    numberOfNodes = len(el)
                    num = math.ceil(math.log(numberOfNodes,1 + self.epsilon/2))
                    num = math.pow(1 + self.epsilon/2,num)
                    i = bisect.bisect_left(GVector, num)
                    if i == len(GVector):
                        i -= 1
                    currentG[i] += 1
                    packer.items.append(Item('A', round(GVector[i])))
        return partition,tuple(currentG)

    def getBestPartitionMethod1(self,listPartitions,listTrees,PossiblePartitionsCost):
        GVector = self.getGVector()
        lastPartitions = [None] * len(listTrees)#np.empty(len(listTrees), dtype=int)
        indexLastPartitions = np.empty(len(listTrees), dtype=int)
        indexLastPartitions.fill(-1)
        for p in listPartitions:
            packer = Binpacker(round((1+self.epsilon)*self.n/self.k))  # Dimension for each bin
            partition,costCurrentPartition = self.createPartitionNodes(listTrees,p,PossiblePartitionsCost,GVector,packer,lastPartitions,indexLastPartitions)
            if(len(packer.items)>=self.k and packer.pack_items(self.k)):       # Return True se riesce a farceli stare
                return partition,p

    def orderFinalPartitions(self,partitions):
        PossiblePartitionsCost = np.array(list(itertools.product(*partitions)))
        PossiblePartitionsCostMerged = PossiblePartitionsCost.sum(axis=1)
        return (np.argsort(PossiblePartitionsCostMerged)),partitions

    def createPartitionNodesMethod2(self,listTrees, i, index, GVector):
        temp = []
        temp.append(listTrees[i])
        currentG = np.zeros(len(GVector))
        tempPartition = self.getPartitionInTree(temp, index, [0])
        if (tempPartition != None):
            for el in tempPartition:
                numberOfNodes = len(el)
                num = math.ceil(math.log(numberOfNodes, 1 + self.epsilon / 2))
                num = math.pow(1 + self.epsilon / 2, num)
                pos = bisect.bisect_left(GVector, num)
                #pos = bisect.bisect_left(GVector, numberOfNodes)
                if pos == len(GVector):
                    pos -= 1
                currentG[pos] += 1
        return currentG  #,coste

    def getBestPartitionMethod2(self,listTrees,possiblePartitionsCost,vectorG,index,GVector,packer,gDictionary):
        if(index == len(listTrees)):
            packer.items = []
            for el in range(len(vectorG)):
                for j in range(int(vectorG[el])):
                    packer.items.append(Item('A', round(GVector[el])))
            vectorG = tuple(vectorG)
            if (vectorG in self.feasibleSets):
                return 0,[]
            elif(vectorG not in self.unfeasibleSets and len(packer.items) >= self.k and packer.pack_items(self.k)):  # Return True se riesce a farceli stare
                self.feasibleSets.add(vectorG)
                return 0,[]
            elif(vectorG not in self.unfeasibleSets):
                self.unfeasibleSets.add(vectorG)
            return math.inf,[]
        costs = []
        partitions = []
        for i in range(len(possiblePartitionsCost[index])):
            currentG = self.createPartitionNodesMethod2(listTrees, index, i, GVector)
            newVectorG = (vectorG+currentG)
            if((index,tuple(newVectorG)) in gDictionary):
                cost, minP = gDictionary[(index,tuple(newVectorG))]
                minP = list(minP)
            else:
                cost, minP = self.getBestPartitionMethod2(listTrees,
                                                        possiblePartitionsCost, (vectorG + currentG),
                                                        index + 1, GVector, packer, gDictionary)
                gDictionary[(index,tuple(newVectorG))] = cost,tuple(minP)
            costs.append(cost+possiblePartitionsCost[index][i])
            partitions.append(minP)
        min = np.argmin(costs)
        partitions[min].append(min)
        return (costs[min],partitions[min])

    def getPartitionAtIPos(self,listTrees,bestP):
        partition = []
        # cost = 0
        for i in range(len(listTrees)):
            temp = []
            temp.append(listTrees[-i - 1])
            tempPartition = self.getPartitionInTree(temp, bestP[i], [0])
            if (tempPartition != None):
                for el in tempPartition:
                    partition.append(el)
        return partition

    def printPartition(self,parition):
        if (parition != None):
            for el in parition:
                print(str(el) + " ")
        else:
            print("non partition exists")

    def getCostPartition(self, parition, G):
        totalCost = 0
        for i in range(len(parition)):
            for j in range(i + 1, len(parition)):
                totalCost += self.calculateCostFinal(G, parition[i], parition[j])
        return totalCost

    def printGraphWithPartitions(self,parition,G):
        pos = nx.spring_layout(G)  # positions for all nodes
        vet = ['r', 'g', 'b', 'y']
        for i in range(len(parition)- 4 ):
            r = lambda: random.randint(0, 255)
            vet.append('#%02X%02X%02X' % (r(), r(), r()))
        for i in range(len(parition)):
            nx.draw_networkx_nodes(G, pos,
                                   nodelist=parition[i],
                                   node_color=vet[i],
                                   node_size=500,
                                   alpha=0.8)
            # edges
        nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
        plt.show()

    def evaluatePartitionsMethod2(self,listTrees,listTotalCosts):
        GVector = self.getGVector()
        vectorG = np.zeros(len(GVector), dtype=int)
        packer = Binpacker(round((1 + self.epsilon) * self.n / self.k))
        gDictionary = dict()
        cost,bestP = self.getBestPartitionMethod2(listTrees, listTotalCosts, vectorG, 0, GVector, packer, gDictionary)
        bestP = self.getPartitionAtIPos(listTrees, bestP)
        return bestP

    def  evaluatePartitionsMethod1(self,listTrees,listTotalCosts):
        fin, PossiblePartitionsCost = self.orderFinalPartitions(
            listTotalCosts)  # get all partition ordered by cost
        bestPBest, p = self.getBestPartitionMethod1(fin, listTrees, listTotalCosts)  # return best partition
        return bestPBest

def saveToFile(file,x,y):
    x = (' '.join(str(e) for e in x))+" time: "
    y = ' '.join(str(e) for e in y)
    file = open(file, "w")
    file.write(x)
    file.write(y)
    file.close()

def testK(samplesForValue,valuesOfK,G,epsilon):
    # G = nx.gaussian_random_partition_graph(50,2,0.1,0.6,0.6)
    # G = nx.connected_watts_strogatz_graph(20,2,0.1)
    x = []
    y = []
    for val in valuesOfK:
        partitor = GraphPartitioning(epsilon, G.number_of_nodes(), val, G)
        for sample in range(samplesForValue):
            if (nx.is_connected(G)):
                root = Tree
                root.data = G
                print("Avvio trasformazione di grafo in albero")
                # partitor.division(root,epsilon*n/(3*k),G)   # create Tree
                start_time = timeit.default_timer()
                partitor.bilancedPartion(root, epsilon * G.number_of_nodes() / (3 * val), G)
                listTrees = partitor.pruneTree(root)
                listTotalCosts = partitor.scanTrees(listTrees)
                # bestP = partitor.evaluatePartitionsMethod1(listTrees, listTotalCosts)
                bestP = partitor.evaluatePartitionsMethod2(listTrees, listTotalCosts)
                elapsed = timeit.default_timer() - start_time
                x.append(val)
                y.append(elapsed)
    print(x)
    print(y)
    saveToFile("test.txt", x, y)

def testEpsilon(samplesForValue,valuesOfEpsilon,G,k):
    # G = nx.gaussian_random_partition_graph(50,2,0.1,0.6,0.6)
    # G = nx.connected_watts_strogatz_graph(20,2,0.1)
    x = []
    y = []
    for val in valuesOfEpsilon:
        partitor = GraphPartitioning(val, G.number_of_nodes(), k, G)
        for sample in range(samplesForValue):
            if (nx.is_connected(G)):
                root = Tree
                root.data = G
                print("Avvio trasformazione di grafo in albero")
                # partitor.division(root,epsilon*n/(3*k),G)   # create Tree
                start_time = timeit.default_timer()
                partitor.bilancedPartion(root, val * G.number_of_nodes() / (3 * k), G)
                listTrees = partitor.pruneTree(root)
                listTotalCosts = partitor.scanTrees(listTrees)
                # bestP = partitor.evaluatePartitionsMethod1(listTrees, listTotalCosts)
                bestP = partitor.evaluatePartitionsMethod2(listTrees, listTotalCosts)
                elapsed = timeit.default_timer() - start_time
                x.append(val)
                y.append(elapsed)
    print(x)
    print(y)
    saveToFile("test.txt", x, y)

def testN(samplesForValue,k,gList,epsilon):
    # G = nx.gaussian_random_partition_graph(50,2,0.1,0.6,0.6)
    # G = nx.connected_watts_strogatz_graph(20,2,0.1)
    x = []
    y = []
    for i in range(len(gList)):
        partitor = GraphPartitioning(epsilon, gList[i].number_of_nodes(), k, gList[i])
        for sample in range(samplesForValue):
            if (nx.is_connected(gList[i])):
                root = Tree
                root.data = gList[i]
                print("Avvio trasformazione di grafo in albero")
                start_time = timeit.default_timer()
                # partitor.division(root,epsilon*n/(3*k),G)   # create Tree
                partitor.bilancedPartion(root, epsilon * gList[i].number_of_nodes() / (3 * k), gList[i])
                listTrees = partitor.pruneTree(root)
                listTotalCosts = partitor.scanTrees(listTrees)
                # bestP = partitor.evaluatePartitionsMethod1(listTrees, listTotalCosts)
                bestP = partitor.evaluatePartitionsMethod2(listTrees, listTotalCosts)
                elapsed = timeit.default_timer() - start_time
                x.append(gList[i].number_of_nodes())
                y.append(elapsed)
    print(x)
    print(y)
    saveToFile("test.txt", x, y)

def main():
    #G = nx.read_weighted_edgelist("graph.csv", delimiter=' ', nodetype=str)
    G = nx.gaussian_random_partition_graph(20,2,0.1,0.6,0.6)
    #G = nx.connected_watts_strogatz_graph(20,2,0.1)
    #G = nx.dorogovtsev_goltsev_mendes_graph(2)
    for (u, v) in G.edges():
        G.edge[u][v]['weight'] = 1#random.randint(0, 10)

    testEpsilon(5,[0.3,0.5],G,3)

    #testK(5, [3, 4], G, 0.5)

    #G2 = nx.gaussian_random_partition_graph(100, 2, 0.1, 0.6, 0.6)
    #for (u, v) in G2.edges():
    #    G2.edge[u][v]['weight'] = 1#random.randint(0, 10)
    #listG = []
    #listG.append(G)
    #listG.append(G2)
    #testN(5,3,[G,G2],0.5)
main()