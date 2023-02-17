import random
import math
import networkx as nx
import igraph as ig
import itertools
import numpy as np
import sys
import bisect

if sys.version_info[0] > 2:
    from src.BinPackerDynamic import Binpacker, Item
    from src.kernighanLin import kernighan_lin_bisection
    from src.Tree import Tree
else:
    from src import *

import matplotlib.pyplot as plt


class GraphPartitioning(object):
    def __init__(self, epsilon, n, k, G):
        self.epsilon = epsilon
        self.n = n
        self.k = k
        self.G = G
        # self.num = 0
        self.root = Tree
        self.root.data = G
        self.num = self.getWeightsOfNodes()

        print("Is it connected?:   " + str(nx.is_connected(G)))
        if nx.is_connected(G) == False:
            print("Removing Singols....")
            listOfNodes = G.nodes()
            for x in range(len(listOfNodes)):
                if G.degree(str(listOfNodes[x])) == 0:
                    G.remove_node(str(listOfNodes[x]))
            print("Is it connected now?:   " + str(nx.is_connected(G)))

    def fromNetworkXtoiGraph(
        self, subTree
    ):  # Conversion from networkX data structure to igraph data
        newlist = []
        edgeList = subTree.data.edges()
        for x in range(
            len(edgeList)
        ):  # Conversion from NetworkX Graph to igraph data structure
            couple = []
            couple.append(
                int(edgeList[x][0])
            )  # the edges has strucure [(n1,n2,{}),....] {} indicates addition attribute
            couple.append(int(edgeList[x][1]))
            newlist.append(couple)
        return ig.Graph(len(subTree.data.nodes()), newlist)

    def createSubTree(self, partition, subTree, threshold, G):  #
        X = G.subgraph(
            str(x) for x in partition[0]
        )  # spotting the two sub-sets divided inside the original Graph G
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
        if (
            len(partition[0]) > threshold
        ):  # progressive recursive approach from root to leaves
            self.division(leftChild, threshold, G)
        if (
            len(partition[1]) > threshold
        ):  # treshold is evaluated in order to avoid the creation of parts that need to be removed  (epsilon*n)/ 3*k
            self.division(rightChild, threshold, G)

    def division(self, subTree, threshold, G):
        if subTree.data.number_of_nodes() <= 1:  # if subtree is a leaf, return it
            return subTree
        else:
            balanceCut = int(
                round(
                    (((self.epsilon / 3) / (1 + self.epsilon / 3)) / 2)
                    * subTree.data.number_of_nodes()
                )
            )  # Evaluation of epsilon' according to formula provided in paper (Section "3. Algortithm and Analysis")
            listOfNodes = (
                subTree.data.nodes()
            )  # Selection of two randon nodes in order to divide the graph in two separate parts
            nodeSource = random.choice(listOfNodes)  # source node
            listOfNodes.remove(nodeSource)
            nodeDestination = random.choice(listOfNodes)  # destination node
            graphNew = self.fromNetworkXtoiGraph(subTree)
            mc = graphNew.mincut(
                int(nodeSource), int(nodeDestination)
            )  # cutting method of igraph library
            partition = mc.partition
            while (
                True
            ):  # igraph has a continus declaration of nodes, so it is not possible to declare a non-continue set of numbers
                cont = 0  # ex: set(3,4) -> creates the set (1,2,3,4) [1,2] needs to be removed
                n = 0
                for x in range(
                    len(partition[0])
                ):  # This "for" fix the problem of igraph continuity (2 line before)
                    if str(partition[0][x - n]) in subTree.data.nodes():
                        cont += 1
                    else:
                        del partition[0][x - n]  # removing the additional elements
                        n += 1
                if (
                    cont >= balanceCut and len(partition[1]) >= balanceCut
                ):  # if the two subsets respect the bipartion constraints of
                    break  # (e')/2*|V| we have created the two subsets.
                else:
                    listOfNodes = (
                        subTree.data.nodes()
                    )  # otherwise we search another bipartion, according the
                    nodeSource = random.choice(listOfNodes)  # constraints
                    listOfNodes.remove(nodeSource)
                    nodeDestination = random.choice(listOfNodes)
                    mc = graphNew.mincut(
                        int(nodeSource), int(nodeDestination)
                    )  # defining a new partition
                    partition = mc.partition
            self.createSubTree(partition, subTree, threshold, G)

    def bilancedPartition(self, subTree, threshold, G):
        if subTree.data.number_of_nodes() <= 1:  # if subtree is a leaf, return it
            return subTree
        else:
            partition = kernighan_lin_bisection(subTree.data)
            X = G.subgraph(
                partition[0]
            )  # spotting the two sub-sets of division inside the original Graph G
            Y = G.subgraph(partition[1])
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
            if (
                self.getWeightsNodes(leftChild.data.node) > threshold
            ):  # comments are the same of createSubTree function (before seen)
                self.bilancedPartition(leftChild, threshold, G)
            if self.getWeightsNodes(rightChild.data.node) > threshold:
                self.bilancedPartition(rightChild, threshold, G)

    def getWeightsNodes(self, nodes):  # calculate of the total weight of nodes
        weight = 0
        for node in nodes:
            weight += nodes[node]["weight"]
        return weight

    def removeBigNodes(self, root, threshold, listTrees):
        if self.getWeightsNodes(root.data.node) < threshold:
            listTrees.append(root)
            return
        else:
            self.removeBigNodes(root.left, threshold, listTrees)
            self.removeBigNodes(root.right, threshold, listTrees)

    def pruneTree(
        self, root
    ):  # remove nodes that are too big by returning the root the the nwe sub trees
        listTrees = []
        num = 0
        for node in self.G.node:
            num += self.G.node[node]["weight"]
        self.removeBigNodes(root, (1 + self.epsilon) * (num / self.k), listTrees)
        return listTrees

    def getGVector(
        self,
    ):  # generate the size for the G vector (paragraph 3.1 of the paper)
        num = 0
        for node in self.G.node:
            num += self.G.node[node]["weight"]
        vector_gValue = []
        val = math.pow(
            (1 + self.epsilon / 2),
            math.floor(
                math.log(self.epsilon * num / (3 * self.k), (1 + self.epsilon / 2))
            ),
        )
        while val < ((1 + self.epsilon) * (num / self.k)):
            vector_gValue.append(val)
            val *= 1 + self.epsilon / 2
        return vector_gValue

    def calculateCost(
        self, graphData, left, right
    ):  # calculate cost of divide two sets
        sum = 0
        leftNodes = left.data.node
        rightNodes = right.data.node
        smallPart = None
        bigPartPart = None
        if len(leftNodes) > len(rightNodes):
            smallPart = rightNodes
            bigPartPart = leftNodes
        else:
            smallPart = leftNodes
            bigPartPart = rightNodes
        for edge in graphData.edge:
            if edge in smallPart:
                edges = graphData.edge[edge]
                for singEdge in edges:
                    if singEdge in bigPartPart:
                        sum += edges[singEdge]["weight"]
        return sum

    def calculateCostFinal(
        self, graphData, left, right
    ):  # calculate final cost of a partitioning
        sum = 0
        leftNodes = left
        rightNodes = right
        smallPart = None
        bigPartPart = None
        if len(leftNodes) > len(rightNodes):
            smallPart = rightNodes
            bigPartPart = leftNodes
        else:
            smallPart = leftNodes
            bigPartPart = rightNodes
        for edge in graphData.edge:
            if edge in smallPart:
                edges = graphData.edge[edge]
                for singEdge in edges:
                    if singEdge in bigPartPart:
                        sum += edges[singEdge]["weight"]  # = 1
        return sum

    def getPartitionInTree(
        self, listNodes, num, i
    ):  # retrieve the num partition in a tree
        tempVett = []
        if i[0] == num:
            temp = []
            for el in listNodes:
                temp.append(el.data.node)
            return temp
        for k in range(len(listNodes)):
            if listNodes[k].left != None:
                father = listNodes[k]
                listNodes.pop(k)
                listNodes.insert(k, father.left)
                listNodes.insert(k + 1, father.right)
                i[0] += 1
                tempVett = self.getPartitionInTree(listNodes, num, i)
                if tempVett != None:
                    return tempVett
                listNodes.pop(k)
                listNodes.pop(k)
                listNodes.insert(k, father)

    def createPossiblePartitions(
        self, listNodes, listTreesCosts, cost
    ):  # generate the cost of each sub-partition in each tree
        listTreesCosts.append(cost)  # append the cost in the cost vector
        for i in range(len(listNodes)):
            if listNodes[i].left != None:
                father = listNodes[i]
                listNodes.pop(i)  # remove the father and add the two sons
                listNodes.insert(i, father.left)
                listNodes.insert(i + 1, father.right)
                newCost = self.calculateCost(father.data, father.left, father.right)
                self.createPossiblePartitions(
                    listNodes, listTreesCosts, cost + newCost
                )  # the cost of the partition is the old cost plus the cost of dividing father.left and father.right
                listNodes.pop(i)
                listNodes.pop(i)
                listNodes.insert(i, father)

    def scanTrees(
        self, listTrees
    ):  # generate the cost of each sub-partition for each sub tree
        listOfTotalCosts = []
        for tree in listTrees:
            temp = []
            listTreesCosts = []
            temp.append(tree)
            self.createPossiblePartitions(temp, listTreesCosts, 0)
            listOfTotalCosts.append(listTreesCosts)
        return listOfTotalCosts

    def getPartitionCost(
        self, n, PossiblePartitionsCost
    ):  # retrieve the cost of the n partition
        cost = 0
        for i in range(len(PossiblePartitionsCost)):
            module = len(PossiblePartitionsCost[-i - 1])
            cost += PossiblePartitionsCost[-i - 1][int(n % module)]
            n -= n % module
            n /= module
        return cost

    def createPartitionNodesMethod1(
        self,
        listTrees,
        n,
        PossiblePartitionsCost,
        GVector,
        packer,
        lastPartitions,
        indexLastPartitions,
    ):  # retrive partition at pos n and return its G vector
        partition = []
        currentG = np.zeros(len(GVector))
        for i in range(len(listTrees)):  # tree in listTrees[::-1]:
            temp = []
            temp.append(listTrees[-i - 1])
            module = len(PossiblePartitionsCost[-i - 1])
            tempPartition = self.getPartitionInTree(temp, n % module, [0])
            indexLastPartitions[i] = n
            lastPartitions[i] = tempPartition
            n -= n % module
            n /= module
            if tempPartition != None:
                for el in tempPartition:
                    partition.append(el)
                    num = math.ceil(math.log(self.num, 1 + self.epsilon / 2))
                    num = math.pow(1 + self.epsilon / 2, num)
                    i = bisect.bisect_left(GVector, num)
                    if i == len(GVector):
                        i -= 1
                    currentG[i] += 1
                    packer.items.append(Item("A", round(GVector[i])))
        return partition, tuple(currentG)

    def getBestPartitionMethod1(
        self, listPartitions, listTrees, PossiblePartitionsCost
    ):
        GVector = self.getGVector()
        lastPartitions = [None] * len(listTrees)  # np.empty(len(listTrees), dtype=int)
        indexLastPartitions = np.empty(len(listTrees), dtype=int)
        indexLastPartitions.fill(-1)
        for p in listPartitions:
            packer = Binpacker(
                round((1 + self.epsilon) * self.n / self.k)
            )  # Dimension for each bin
            partition, costCurrentPartition = self.createPartitionNodesMethod1(
                listTrees,
                p,
                PossiblePartitionsCost,
                GVector,
                packer,
                lastPartitions,
                indexLastPartitions,
            )
            if len(packer.items) >= self.k and packer.pack_items(
                self.k
            ):  # Return True if the binPackingProblem is satisfied
                return partition, p

    def orderFinalPartitions(self, partitions):
        PossiblePartitionsCost = np.array(list(itertools.product(*partitions)))
        PossiblePartitionsCostMerged = PossiblePartitionsCost.sum(axis=1)
        return (np.argsort(PossiblePartitionsCostMerged)), partitions

    def createPartitionNodesMethod2(
        self, listTrees, i, index, GVector
    ):  # generate the G vector for the element
        temp = []
        temp.append(listTrees[i])
        currentG = np.zeros(len(GVector))
        tempPartition = self.getPartitionInTree(temp, index, [0])
        if tempPartition != None:
            for el in tempPartition:
                weightOfNodes = self.getWeightsNodes(el)
                num = math.ceil(math.log(weightOfNodes, 1 + self.epsilon / 2))
                num = math.pow(1 + self.epsilon / 2, num)
                pos = bisect.bisect_left(GVector, num)
                # pos = bisect.bisect_left(GVector, numberOfNodes)
                if pos == len(GVector):
                    pos -= 1
                currentG[pos] += 1
        return currentG  # ,coste

    def getBestPartitionMethod2(
        self,
        listTrees,
        possiblePartitionsCost,
        vectorG,
        index,
        GVector,
        packer,
        gDictionary,
    ):  # dynamic programming part
        if index == len(listTrees):  # if I gerate a partition verify if it is feasible
            packer.items = []
            for el in range(len(vectorG)):
                for j in range(int(vectorG[el])):
                    packer.items.append(Item("A", round(GVector[el])))
            if len(packer.items) >= self.k and packer.pack_items(
                self.k
            ):  # BinPackingProblem(we use the BinPackerDynamic library)
                return 0, []  # return cost = 0 if feasible
            return float("inf"), []  # return cost = infinite if unfeasible
        costs = []
        partitions = []
        for i in range(len(possiblePartitionsCost[index])):
            currentG = self.createPartitionNodesMethod2(
                listTrees, index, i, GVector
            )  # create the G vector for the selected element
            newVectorG = (
                vectorG + currentG
            )  # cumulate the new G vetor with the one generated in the previous tree
            if (
                index,
                tuple(newVectorG),
            ) in gDictionary:  # if I already know the result, stop the recursion and return the result
                cost, minP = gDictionary[(index, tuple(newVectorG))]
                minP = list(minP)
            else:  # if I still don' t know the result, calculate it
                cost, minP = self.getBestPartitionMethod2(
                    listTrees,
                    possiblePartitionsCost,
                    (vectorG + currentG),
                    index + 1,
                    GVector,
                    packer,
                    gDictionary,
                )
                gDictionary[(index, tuple(newVectorG))] = cost, tuple(
                    minP
                )  # save the result in the dictionary for the future
            costs.append(cost + possiblePartitionsCost[index][i])
            partitions.append(minP)
        min = np.argmin(costs)
        partitions[min].append(min)
        return (
            costs[min],
            partitions[min],
        )  # return the temp-partition with the lowest cost

    def getPartitionAtIPos(
        self, listTrees, bestP
    ):  # convert number of partition in a set of elements
        partition = []
        for i in range(len(listTrees)):
            temp = []
            temp.append(listTrees[-i - 1])
            tempPartition = self.getPartitionInTree(
                temp, bestP[i], [0]
            )  # get partition for each tree
            if tempPartition != None:
                for el in tempPartition:
                    partition.append(el)
        return partition

    def printPartition(self, parition):  # print textually the best partition
        if parition != None:
            for el in parition:
                print(str(el) + " ")
        else:
            print("no partitions exist")

    def getCostPartition(
        self, parition
    ):  # calculate and print the cost of the best partition
        totalCost = 0
        for i in range(
            len(parition)
        ):  # for each couple of subsets in the partition calculate the cost
            for j in range(i + 1, len(parition)):
                totalCost += self.calculateCostFinal(self.G, parition[i], parition[j])
        return totalCost  # return the cost for divide all the subsets in the best partition

    def printGraphWithPartitions(
        self, parition
    ):  # render graphically the best partitions
        pos = nx.spring_layout(self.G)  # positions for all nodes
        vet = ["r", "g", "b", "y"]  # array with all the colors
        for i in range(len(parition) - len(vet)):
            r = lambda: random.randint(0, 255)  # generate random color
            vet.append("#%02X%02X%02X" % (r(), r(), r()))
        for i in range(len(parition)):  # render each partition with the desired color
            nx.draw_networkx_nodes(
                self.G,
                pos,
                nodelist=parition[i],
                node_color=vet[i],
                node_size=500,
                alpha=0.8,
            )
            # edges
        nx.draw_networkx_edges(self.G, pos, width=1.0, alpha=0.5)  # render the edges
        plt.show()

    def getWeightsOfNodes(self):  # calculate sum of all the weights of the nodes
        num = 0
        for node in self.G.node:
            num += self.G.node[node]["weight"]
        self.num = num
        return num

    def evaluatePartitionsMethod2(
        self, listTrees, listTotalCosts
    ):  # method that uses dynamic programming (paragraph 3.2 of the paper)
        GVector = self.getGVector()
        vectorG = np.zeros(len(GVector), dtype=int)
        packer = Binpacker(round((1 + self.epsilon) * self.num / self.k))
        gDictionary = dict()
        cost, bestP = self.getBestPartitionMethod2(
            listTrees, listTotalCosts, vectorG, 0, GVector, packer, gDictionary
        )
        bestP = self.getPartitionAtIPos(
            listTrees, bestP
        )  # convert number of partition in a set of elements
        return bestP

    def evaluatePartitionsMethod1(
        self, listTrees, listTotalCosts
    ):  # method that doesn' t use dynamic programming
        fin, PossiblePartitionsCost = self.orderFinalPartitions(
            listTotalCosts
        )  # generate each partition ordered by cost
        bestPBest, p = self.getBestPartitionMethod1(
            fin, listTrees, listTotalCosts
        )  # return best partition (lowest cost that is feasible)
        return bestPBest

    def run(self):
        # self.division(self.root, self.epsilon * self.num / (3 * self.k), self.G)
        self.bilancedPartition(
            self.root, self.epsilon * self.num / (3 * self.k), self.G
        )  # partition the graph and generate the tree
        listTrees = self.pruneTree(self.root)  # prune the tree
        listTotalCosts = self.scanTrees(
            listTrees
        )  # calculate the cost for each partition inside each tree

        # select the desired method of execution
        # bestP = partitor.evaluatePartitionsMethod1(listTrees, listTotalCosts)    # method that doesn' t use dynamic programming
        bestP = self.evaluatePartitionsMethod2(
            listTrees, listTotalCosts
        )  # method that uses dynamic programming

        return bestP

    def savePartitionToFile(self, bestP, file):
        stringToSave = ""
        for el in bestP:
            stringToSave += "["
            for key in el.keys():
                stringToSave += str(key) + ","
            stringToSave = stringToSave[:-1]
            stringToSave += "]"
        file = open("bestPartition.txt", "w")
        file.write(stringToSave)
        file.close()


def main(sourceFile, destinationFile, k, epsilon):
    if sourceFile != "":
        G = nx.read_weighted_edgelist(
            sourceFile, delimiter=" ", nodetype=str
        )  # read the graph from file(adjacency list) graph.csv is an example of format
    else:
        G = nx.gaussian_random_partition_graph(
            100, 2, 0.1, 0.6, 0.6
        )  # select and generate the intial graph
        # G = nx.connected_watts_strogatz_graph(15,2,0.1)
        # G = nx.dorogovtsev_goltsev_mendes_graph(3)
        for u, v in G.edges():  # initialize the weight of the edges of the graph
            G.edge[u][v]["weight"] = 1  # random.randint(0, 10)
    for n in G.nodes():  # initialize the weight of the nodes of the graph
        G.node[n]["weight"] = 1  # random.random() * 100
    n = G.number_of_nodes()
    if nx.is_connected(G) and k <= n:
        print("Avvio trasformazione di grafo in albero")
        partitor = GraphPartitioning(
            epsilon, n, k, G
        )  # create the object and pass it the intial graph
        bestP = partitor.run()  # run the algorithm and return the best paritition
        print("The best parition is:")
        partitor.printPartition(bestP)  # print textually the best partition
        print("The cost of the parition is:")
        print(
            partitor.getCostPartition(bestP)
        )  # calculate and print the cost of the best partition
        partitor.printGraphWithPartitions(
            bestP
        )  # render graphically the best partitions
        partitor.savePartitionToFile(
            bestP, destinationFile
        )  # save the best partition to file


if len(sys.argv) == 1:
    main(
        "", "bestPartition.txt", 3, 0.9
    )  # if sourceFile=="" a networkx graph is generated otherwise the file is read from graph
    # main("graph.csv","bestPartition.txt",3,0.9)
else:
    main(
        sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
    )  # read input from command line parameters
