import sys

if sys.version_info[0] == 2:
    from BinPackerDynamic import Binpacker, Item
    from kernighanLin import kernighan_lin_bisection
    from src.Tree import Tree, TreeUtil
