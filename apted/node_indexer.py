#
# The MIT License
#
# Copyright 2017 Joao Felipe Pimentel

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
"""NodeIndexer and NodeInfo implementations"""
from __future__ import (absolute_import, division)

from .helpers import Config

class NodeIndexer(object):
    """Indexes nodes of the input tree to the algorithm that is already
    parsed to tree structure. Stores various indices on nodes required for
    efficient computation of APTED [1,2]. Additionally, it stores single-value
    properties of the tree.

    For indexing we use four tree traversals that assign ids to the nodes:

    - left-to-right preorder [1],
    - right-to-left preorder [1],
    - left-to-right postorder [2],
    - right-to-left postorder [2].

    See the source code for more algorithm-related comments.

    References:
    [1] M. Pawlik and N. Augsten. Efficient Computation of the Tree Edit
        Distance. ACM Transactions on Database Systems (TODS) 40(1). 2015.
    [2] M. Pawlik and N. Augsten. Tree edit distance: Robust and memory-
        efficient. Information Systems 56. 2016.
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self, tree, config=None):
        # pylint: disable=too-many-statements

        self.config = config or Config()
        """Config object that specifies how to calculate the edit distance"""

        self.preorder_tmp = 0
        """Temporary variable that stores preorder index"""

        self.pre_ltr_info = []
        """Map left-to-right preorder index to NodeInfo"""

        self.post_ltr_info = []
        """Map left-to-right postorder index to NodeInfo"""

        self.post_rtl_info = []
        """Map right-to-left postorder index to NodeInfo"""

        self.pre_rtl_info = []
        """Map right-to-left preorder index to NodeInfo"""

        self.tree_size = 0
        """Size of input tree"""

        self.lchl = 0
        """Number of leftmost-child leaf nodes in the input tree
        See [2, Section 5.3]."""

        self.rchl = 0
        """Number of rightmost-child leaf nodes in the input tree
        See [2, Section 5.3]."""

        self.current_index = 0
        """Left-to-right preorder id of the current subtree's root node.
        Used in the tree decomposition phase of APTED. See [1, Algorithm 1]"""

        root, _ = self.index_nodes(tree, -1)
        root.parent = NodeInfo(None, -1)
        root.parent.fake_child = root


        self.tree_size = self.pre_ltr_info[0].size
        self.pre_rtl_info = [None] * self.tree_size
        self.post_rtl_info = [None] * self.tree_size

        self.post_traversal_indexing()

        # Fallback
        self.parents = []
        self.children = []
        self.node_type_l = []
        self.node_type_r = []
        self.pre_ltr_to_desc_sum = []
        self.pre_ltr_to_kr_sum = []
        self.pre_ltr_to_rev_kr_sum = []
        self.pre_ltr_to_node = []
        self.sizes = []
        self.pre_ltr_to_pre_rtl = []
        self.pre_ltr_to_post_ltr = []
        self.pre_ltr_to_post_rtl = []
        self.pre_ltr_to_ln = []
        self.pre_ltr_to_sum_del_cost = []
        self.pre_ltr_to_sum_ins_cost = []

        for node_info in self.pre_ltr_info:
            self.parents.append(getattr(node_info.parent, 'pre_ltr', -1))
            self.children.append([x.pre_ltr for x in node_info.children])
            self.node_type_l.append(node_info.type_l)
            self.node_type_r.append(node_info.type_r)
            self.pre_ltr_to_desc_sum.append(node_info.desc_sum)
            self.pre_ltr_to_kr_sum.append(node_info.kr_sum)
            self.pre_ltr_to_rev_kr_sum.append(node_info.rev_kr_sum)
            self.pre_ltr_to_node.append(node_info.node)
            self.sizes.append(node_info.size)
            self.pre_ltr_to_pre_rtl.append(node_info.pre_rtl)
            self.pre_ltr_to_post_ltr.append(node_info.post_ltr)
            self.pre_ltr_to_post_rtl.append(node_info.post_rtl)
            self.pre_ltr_to_ln.append(getattr(node_info.lnl, 'pre_ltr', -1))
            self.pre_ltr_to_sum_del_cost.append(node_info.sum_del_cost)
            self.pre_ltr_to_sum_ins_cost.append(node_info.sum_ins_cost)
        self.pre_rtl_to_pre_ltr = []
        self.pre_rtl_to_ln = []
        for node_info in self.pre_rtl_info:
            self.pre_rtl_to_pre_ltr.append(node_info.pre_ltr)
            self.pre_rtl_to_ln.append(getattr(node_info.lnr, 'pre_rtl', -1))
        self.post_ltr_to_pre_ltr = []
        self.post_ltr_to_lld = []
        for node_info in self.post_ltr_info:
            self.post_ltr_to_pre_ltr.append(node_info.pre_ltr)
            self.post_ltr_to_lld.append(node_info.lld.post_ltr)
        self.post_rtl_to_pre_ltr = []
        self.post_rtl_to_rld = []
        for node_info in self.post_rtl_info:
            self.post_rtl_to_pre_ltr.append(node_info.pre_ltr)
            self.post_rtl_to_rld.append(node_info.rld.post_rtl)

        """
        print("Indexes")
        print("P", self.parents)
        print("C", self.children)
        print("NTL", self.node_type_l)
        print("NTR", self.node_type_r)
        print("LDS", self.pre_ltr_to_desc_sum)
        print("LKR", self.pre_ltr_to_kr_sum)
        print("LRKR", self.pre_ltr_to_rev_kr_sum)
        print("LN", self.pre_ltr_to_node)
        print("S", self.sizes)
        print("LR", self.pre_ltr_to_pre_rtl)
        print("RL", self.pre_rtl_to_pre_ltr)
        print("PLL", self.post_ltr_to_pre_ltr)
        print("LPL", self.pre_ltr_to_post_ltr)
        print("LPR", self.pre_ltr_to_post_rtl)
        print("PRL", self.post_rtl_to_pre_ltr)
        print("-")
        print("PLlld", self.post_ltr_to_lld)
        print("PRrld", self.post_rtl_to_rld)
        print("LLN", self.pre_ltr_to_ln)
        print("RLN", self.pre_rtl_to_ln)
        print("LSDC", self.pre_ltr_to_sum_del_cost)
        print("LSIC", self.pre_ltr_to_sum_ins_cost)
        print(self.lchl, self.rchl)
        """

    def index_nodes(self, node, postorder):
        """Preprocesses each node of the tree and creates associated NodeInfo
        Computes the following attributes: node, pre_ltr, post_ltr, parent,
          children, size, desc_sum, kr_sum, rev_kr_sum
        Also computes the following indices: pre_ltr_info, post_ltr_info

        It is a recursive method that traverses the tree once

        node is the current node while traversing the input tree
        postorder is the postorder id of the current node

        return NodeInfo for node and accumulated sum of subtree sizes
        rooted at descendant nodes
        """
        desc_sizes = current_size = kr_sizes_sum = revkr_sizes_sum = 0

        node_info = NodeInfo(node, self.preorder_tmp)
        self.preorder_tmp += 1
        self.pre_ltr_info.append(node_info)

        children = self.config.children(node)
        last_index = len(children) - 1

        for child_index, child in enumerate(children):
            child_info, desc_sizes_tmp = self.index_nodes(child, postorder)
            child_info.parent = node_info
            node_info.children.append(child_info)

            postorder = child_info.post_ltr

            current_size += child_info.size
            desc_sizes += desc_sizes_tmp

            kr_sizes_sum += child_info.kr_sum
            revkr_sizes_sum += child_info.rev_kr_sum
            if child_index == 0:
                kr_sizes_sum -= child_info.size
                child_info.type_l = True
            if child_index == last_index:
                revkr_sizes_sum -= child_info.size
                child_info.type_r = True

        postorder += 1
        self.post_ltr_info.append(node_info)
        node_info.post_ltr = postorder

        node_info.size = nsize = current_size + 1
        desc_sizes_tmp = desc_sizes + nsize

        node_info.desc_sum = (nsize * (nsize + 3)) / 2 - desc_sizes_tmp
        node_info.kr_sum = kr_sizes_sum + nsize
        node_info.rev_kr_sum = revkr_sizes_sum + nsize

        return node_info, desc_sizes_tmp

    def post_traversal_indexing(self):
        """Indexes the nodes of the input tree.
        It computes the following indices, which could not be computed
        immediately while traversing the tree in index_nodes: pre_ltr_to_ln,
          post_ltr_to_lld, post_rtl_to_rlf, pre_rtl_to_ln

        Runs in linear time in the input tree size.
        Currently requires two loops over input tree nodes.
        Can be reduced to one loop (See the code)
        """

        current_leaf = NodeInfo.EMPTY
        node_for_sum = -1
        parent_for_sum = -1

        delete = self.config.delete
        insert = self.config.insert

        for i, node in enumerate(self.pre_ltr_info):
            node.pre_rtl = self.tree_size - 1 - node.post_ltr
            node.post_rtl = self.tree_size - 1 - node.pre_ltr
            self.pre_rtl_info[node.pre_rtl] = node
            self.post_rtl_info[node.post_rtl] = node

            node.lnl = current_leaf
            if not node.children:
                current_leaf = node

                # Count lchl and rchl if node is leaf
                # Note that it must visit parent before child in order
                # to have pre_rtl computed
                parent = node.parent
                if parent:
                    if parent.pre_ltr + 1 == node.pre_ltr:
                        self.lchl += 1
                    elif parent.pre_rtl + 1 == node.pre_rtl:
                        self.rchl += 1

            # Sum up costs of deleting and inserting entire subtrees.
            # Reverse the node index.
            # Traverses nodes bottom-up
            node_for_sum = self.pre_ltr_info[self.tree_size - i - 1]
            parent_for_sum = node_for_sum.parent
            # Update myself
            node_for_sum.sum_del_cost += delete(node_for_sum.node)
            node_for_sum.sum_ins_cost += insert(node_for_sum.node)

            if parent_for_sum:
                parent_for_sum.sum_del_cost += node_for_sum.sum_del_cost
                parent_for_sum.sum_ins_cost += node_for_sum.sum_ins_cost


        current_leaf = NodeInfo.EMPTY
        for i in range(self.tree_size):
            # right-to-left preorder traversal
            node = self.pre_rtl_info[i]
            node.lnr = current_leaf
            if not node.children:
                current_leaf = node

            # left-to-right postorder traversal
            # Stores leftmost leaf descendants for each node.
            # Used for mapping computation.
            node = self.post_ltr_info[i]
            node.lld = node if node.size == 1 else node.children[0].lld

            # right-to-left postorder traversal
            # Stores rightmost leaf descendants for each node.
            node = self.post_rtl_info[i]
            node.rld = node if node.size == 1 else node.children[-1].rld

    def preorder_ltr(self, tree, target=None):
        """Generator that traverses tree in left-to-right preorder"""
        if target is None:
            target = tree.pre_ltr + tree.size
        for pre_ltr in range(tree.pre_ltr, target):
            yield pre_ltr, self.pre_ltr_info[pre_ltr]

    def preorder_rtl(self, tree):
        """Generator that traverses tree in right-to-left preorder"""
        for pre_rtl in range(tree.pre_rtl + tree.size - 1, tree.pre_rtl - 1, -1):
            yield pre_rtl, self.pre_rtl_info[pre_rtl]

    def traverse_up(self, current, target=None):
        """Traverse up to parent until it reaches the target
        If target is not specified, it is the root of the tree
        Generates:
           parent Info
           last visited node
        """
        target = target or self.pre_ltr_info[0]
        target_pre_ltr = target.pre_ltr
        parent = current.parent
        while parent and parent.pre_ltr >= target_pre_ltr:
            yield parent, current
            current, parent = parent, parent.parent

    def walk_up(self, current, target=None):
        """Same thing as traverse_up, but has a initial iteration
        that generates current, None"""
        child = NodeInfo(None, -1)
        child.parent = current
        current.fake_child = child
        yield current, child
        current.fake_child = None
        for parent, curr in self.traverse_up(current, target):
            yield parent, curr

    @property
    def current_node(self):
        return self.pre_ltr_info[self.current_index]


class NodeInfo(object):
    """Represents a Tree Node with extra information"""
    # pylint: disable=too-many-instance-attributes

    EMPTY = None

    def __init__(self, node, preorder):
        if node:
            node.index = preorder

        self.node = node
        """Node referred by this info"""

        self.pre_ltr = preorder
        """Left-to-right preorder traversal index"""

        self.pre_rtl = -1
        """Right-to-left preorder traversal index"""

        self.post_ltr = -1
        """Left-to-right postorder traversal index"""

        self.post_rtl = -1
        """Rigth-to-left postorder traversal index"""

        self.parent = None
        """Parent node_info"""

        self._children = []
        """Node children in left-to-right preorder"""

        self.type_l = False
        """Node lies on the leftmost path starting at its parent
        See [2, Section 5.3, Algorithm 1, Lines 26,36]"""

        self.type_r = False
        """Node lies on the rightmost path starting at its parent
        See [2, Section 5.3, Algorithm 1, Lines 26,36]"""

        self.desc_sum = 0
        """Cost of spf_A (single path function using an inner path)
        for the subtree rooted at this node
        See [1, Section 5.2]
        """

        self.kr_sum = 0
        """Cost of spf_L (single path function using the leftmost path)
        for the subtree rooted at this node
        See [1, Section 5.2]"""

        self.rev_kr_sum = 0
        """Cost of spf_R (single path function using the rightmost path)
        for the subtree rooted at this node
        See [1, Section 5.2]"""

        self.size = 1
        """Size of node subtree (including itself and all its descendants)"""

        self.lnl = self.EMPTY
        """First leaf node to the left of this node
        See [1, Section 8.4]."""

        self.lnr = self.EMPTY
        """First leaf node to the right of this node
        See [1, Section 8.4]."""

        self.lnc = self.EMPTY
        self.fnc = self.EMPTY
        self.ftc = self.EMPTY
        self.ln_in_use = None
        """Current ln node. It changes during the execution"""


        self.lld = self.EMPTY
        """Leftmost leaf descendant of this node"""

        self.rld = self.EMPTY
        """Rightmost leaf descendant of this node"""

        self.sum_del_cost = 0
        """Cost of deleting all nodes in the subtree rooted at this node"""

        self.sum_ins_cost = 0
        """Cost of inserting all nodes in the subtree rooted at this node"""

        self.fake_child = None
        """Overrides leftmost and rightmost child"""


    def __bool__(self):
        return bool(self.node)

    @property
    def children(self):
        """Returns node children"""
        if self.fake_child is not None:
            return [self.fake_child]
        return self._children

    @property
    def rightmost(self):
        """Returns rightmost child"""
        children = self.children
        if not children:
            return NodeInfo(None, -1)
        return children[-1]

    @property
    def leftmost(self):
        """Returns leftmost child"""
        children = self.children
        if not children:
            return NodeInfo(None, -1)
        return children[0]

    def __repr__(self):
        return str(self.pre_ltr)

EMPTY = NodeInfo.EMPTY = NodeInfo(None, -1)
EMPTY.lnl = EMPTY.lnr = EMPTY.lld = EMPTY.rld = EMPTY
