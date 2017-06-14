class Tree(object):                                                         # Tree data structure
    def __init__(self):
        self.left = None
        self.right = None
        self.data = None

class TreeUtil(object):
    @staticmethod
    def getDepth(root):
        if (root==None):
            return 0
        left = TreeUtil.getDepth(root.left)
        right = TreeUtil.getDepth(root.right)
        if(left>right):
            return (left+1)
        return (right+1)

    @staticmethod
    def getNumberOfNodes(root):
        if (root==None):
            return 0
        #left = getNumberOfNodes(root.left)
        #right = getNumberOfNodes(root.right)
        #return left+right