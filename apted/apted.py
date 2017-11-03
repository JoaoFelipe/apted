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
"""APTED implementation"""
from __future__ import (absolute_import, division)

import math
from numpy import zeros
from .helpers import Config
from .node_indexer import NodeIndexer
from .single_path_functions import spf1, SinglePathFunction, LEFT, RIGHT, INNER
# pylint: disable=invalid-name

class Cost(object):
    """Represents a Cost for opt strategy calculation"""
    # pylint: disable=too-few-public-methods
    def __init__(self):
        self.l = 0.0
        self.r = 0.0
        self.i = 0.0
        self.path = 0

    def set_lri(self, value):
        """Sets l, r, and i values"""
        self.l = self.r = self.i = value


class APTED(object):
    """Implements APTED algorithm [1,2].

    - Optimal strategy with all paths.
    - Single-node single path function supports currently only unit cost.
    - Two-node single path function not included.
    - \\Delta^L and \\Delta^R based on Zhang and Shasha's algorithm for executing
      left and right paths (as in [3]). If only left and right paths are used
      in the strategy, the memory usage is reduced by one quadratic array.
    - For any other path \\Delta^A from [1] is used.

    References:
    - [1] M. Pawlik and N. Augsten. Efficient Computation of the Tree Edit
      Distance. ACM Transactions on Database Systems (TODS) 40(1). 2015.
    - [2] M. Pawlik and N. Augsten. Tree edit distance: Robust and memory-
      efficient. Information Systems 56. 2016.
    - [3] M. Pawlik and N. Augsten. RTED: A Robust Algorithm for the Tree Edit
      Distance. PVLDB 5(4). 2011.
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self, tree1, tree2, config=None):
        self.config = config or Config()
        """Config object that specifies how to calculate the edit distance"""

        self.delta = []
        """The distance matrix [1, Sections 3.4,8.2,8.3]
        Used to store intermediate distances between pairs of subtrees"""

        self.q = []
        """One of distance arrays to store intermediate distances in spfA.
        TODO: Verify if other spf-local arrays are initialised within spf.
        If yes, move q to spf to - then, an offset has to be used to access it
        """

        self.fn = []
        """Array used in the algorithm before [1]. Using it does not change
        the complexity.
        TODO: Do not use it [1, Section 8.4]."""

        self.ft = []
        """Array used in the algorithm before [1]. Using it does not change
        the complexity.
        TODO: Do not use it [1, Section 8.4]."""

        self.counter = 0
        """Stores the number of subproblems encountered while computing the
        distance
        See [1, Section 10]."""

        self.it1 = NodeIndexer(tree1, self.config)
        """Stores the indexes of the first input tree"""

        self.it2 = NodeIndexer(tree2, self.config)
        """Stores the indexes of the second input tree"""

        self.result = None


    def compute_edit_distance(self):
        """Compute tree edit distance between source and destination trees
        using APTED algorithm [1,2]."""
        # Initialize delta array
        if self.result is None:
            if self.it1.lchl < self.it1.rchl:
                self.delta = self.compute_opt_strategy_post_l()
            else:
                self.delta = self.compute_opt_strategy_post_r()

            self.ted_init()
            self.result = self.gted()
        return self.result

    def compute_edit_distance_spf_test(self, spf_type):
        """This method is only for testing purspose. It computes TED with a
        fixed path type in the strategy to trigger execution of a specific
        single-path function"""
        # Initialise delta array.
        if self.result is None:
            index_1 = self.it1.pre_ltr_info
            size1, size2 = self.it1.tree_size, self.it2.tree_size
            self.delta = zeros((size1, size2), float)
            # Fix a path type to trigger specific spf.
            for i in range(size1):
                for j in range(size2):
                    # Fix path type
                    if spf_type == LEFT:
                        self.delta[i][j] = index_1[i].lld.pre_ltr + 1
                    elif spf_type == RIGHT:
                        self.delta[i][j] = index_1[i].rld.pre_ltr + 1
            self.ted_init()
            self.result = self.gted()
        return self.result

    def ted_init(self):
        """After the optimal strategy is computed, initializes distances of
        deleting and inserting subtrees without their root nodes."""
        it1, it2 = self.it1, self.it2

        delta = self.delta
        # Reset the subproblems counter.
        self.counter = 0
        # Initialize arrays
        max_size = max(it1.tree_size, it2.tree_size) + 1
        # todo: Move q initialzation to spfA
        self.q = [0.0] * max_size
        # todo: Do not use fn and ft arrays [1, Section 8.4]
        self.fn = [0] * (max_size + 1)
        self.ft = [0] * (max_size + 1)
        # Computer subtree distances without the root nodes when one of
        # the subtrees is a single node

        insert = self.config.insert
        delete = self.config.delete

        # Loop over the nodes in reversed(?) left-to-right preorder.
        for x, node1 in enumerate(it1.pre_ltr_info):
            size1 = node1.size
            for y, node2 in enumerate(it2.pre_ltr_info):
                # Set values in delta based on the sums of deletion and
                # insertion costs. Substract the costs for root nodes.
                # In this method we don't have to verify the order of the
                # input trees because it is equal to the original.
                size2 = node2.size
                if size1 == 1 and size2 == 1:
                    delta[x][y] = 0.0
                elif size1 == 1:
                    delta[x][y] = node2.sum_ins_cost - insert(node2.node)
                elif size2 == 1:
                    delta[x][y] = node1.sum_del_cost - delete(node1.node)

    def compute_opt_strategy_post_l(self):
        """Compute the optimal strategy using left-to-right postorder traversal
        of the nodes. See [2, Algorithm 1].
        """
        def update_parent(min_cost, node, cost, parent_cost):
            """Update parent cost according to node cost and min_cost"""
            update_path = None
            parent_cost.r += min_cost
            tmp_cost = -min_cost + cost.i
            if tmp_cost < parent_cost.i:
                parent_cost.i = tmp_cost
                update_path = parent_cost.path = cost.path
            if node.type_r:
                parent_cost.i += parent_cost.r
                parent_cost.r += cost.r - min_cost
            if node.type_l:
                parent_cost.l += cost.l
            else:
                parent_cost.l += min_cost

            return update_path
        order1 = self.it1.post_ltr_info
        order2 = self.it2.post_ltr_info
        cost_index = lambda node: node.post_ltr
        return self.compute_opt_strategy_post(
            order1, order2, cost_index, update_parent
        )

    def compute_opt_strategy_post_r(self):
        """Compute the optimal strategy using right-to-left postorder traversal
        of the nodes. See [2, Algorithm 1].
        """
        def update_parent(min_cost, node, cost, parent_cost):
            """Update parent cost according to node cost and min_cost"""
            update_path = None
            parent_cost.l += min_cost
            tmp_cost = -min_cost + cost.i
            if tmp_cost < parent_cost.i:
                parent_cost.i = tmp_cost
                update_path = parent_cost.path = cost.path
            if node.type_l:
                parent_cost.i += parent_cost.l
                parent_cost.l += cost.l - min_cost
            if node.type_r:
                parent_cost.r += cost.r
            else:
                parent_cost.r += min_cost

            return update_path
        order1 = self.it1.post_rtl_info
        order2 = self.it2.post_rtl_info
        cost_index = lambda node: node.post_rtl
        return self.compute_opt_strategy_post(
            order1, order2, cost_index, update_parent
        )

    def compute_opt_strategy_post(self, order1, order2, costi, update_parent):
        """Compute the optimal strategy using generic post-ortder traversal
        See [2, Algorithm 1].
        """
        # pylint: disable=too-many-statements
        # pylint: disable=too-many-branches
        # pylint: disable=no-self-use
        # pylint: disable=too-many-locals
        # pylint: disable=cell-var-from-loop
        it1, it2 = self.it1, self.it2
        size1, size2 = it1.tree_size, it2.tree_size
        strategy = zeros((size1, size2), float)
        cost1 = [None] * size1

        leaf_row = [Cost() for _ in range(size2)]

        path_id_offset = size1
        min_cost = float('inf')
        strategy_path = -1

        rows_to_reuse = []
        pre_rtl_1 = it1.pre_rtl_info
        pre_rtl_2 = it2.pre_rtl_info

        for node1 in order1:
            v_cost = costi(node1)
            v_in_pre_ltr = node1.pre_ltr

            strategy_v = strategy[v_in_pre_ltr]

            parent1 = node1.parent
            size_v = node1.size
            kr_sum_v = node1.kr_sum
            rev_kr_sum_v = node1.rev_kr_sum
            desc_sum_v = node1.desc_sum

            # this is the left path's ID which is the leftmost leaf node:
            # l-r_preorder(r-l_preorder(v) + |Fv| - 1)
            left_path_v = -(pre_rtl_1[node1.pre_rtl + size_v - 1].pre_ltr + 1)
            # this is the right path's ID which is the rightmost leaf node:
            # l-r_preorder(v) + |Fv| - 1
            right_path_v = v_in_pre_ltr + size_v


            if not node1.children:
                cost_pointer_v = cost1[v_cost] = leaf_row
                for node2 in order2:
                    w_cost, w_pre = costi(node2), node2.pre_ltr
                    strategy_v[w_pre] = cost_pointer_v[w_cost].path = v_in_pre_ltr
            else:
                cost_pointer_v = cost1[v_cost]

            if parent1:
                parent_v = costi(parent1)
                if cost1[parent_v] is None:
                    if rows_to_reuse:
                        cost1[parent_v] = rows_to_reuse.pop()
                    else:
                        cost1[parent_v] = [Cost() for _ in range(size2)]

                cost_pointer_parent_v = cost1[parent_v]

            cost2 = [Cost() for _ in range(size2)]

            for node2 in order2:
                w_cost = costi(node2)
                w_in_pre_ltr = node2.pre_ltr

                cost_pointer_w = cost2[w_cost]
                cost_pointer_vw = cost_pointer_v[w_cost]

                parent2 = node2.parent

                size_w = node2.size
                if not node2.children:
                    cost_pointer_w.set_lri(0)
                    cost_pointer_w.path = w_in_pre_ltr

                if size_v <= 1 or size_w <= 1:
                    min_cost = max(size_v, size_w)
                    strategy_path = -1
                else:
                    min_cost, _, strategy_fn = min(
                        (
                            size_v * node2.kr_sum + cost_pointer_vw.l, 1,
                            lambda: left_path_v
                        ),
                        (
                            size_v * node2.rev_kr_sum + cost_pointer_vw.r, 2,
                            lambda: right_path_v
                        ),
                        (
                            size_v * node2.desc_sum + cost_pointer_vw.i, 3,
                            lambda: cost_pointer_vw.path + 1
                        ),
                        (
                            size_w * kr_sum_v + cost_pointer_w.l, 4,
                            lambda: -(pre_rtl_2[node2.pre_rtl + size_w - 1].pre_ltr +
                                      path_id_offset + 1)
                        ),
                        (
                            size_w * rev_kr_sum_v + cost_pointer_w.r, 5,
                            lambda: w_in_pre_ltr + size_w + path_id_offset
                        ),
                        (
                            size_w * desc_sum_v + cost_pointer_w.i, 6,
                            lambda: cost_pointer_w.path + path_id_offset + 1
                        )
                    )
                    strategy_path = strategy_fn()

                if parent1:
                    new_path = update_parent(
                        min_cost, node1, cost_pointer_vw,
                        cost_pointer_parent_v[w_cost]
                    )
                    if new_path is not None:
                        strategy_v[w_in_pre_ltr] = new_path

                if parent2:
                    update_parent(
                        min_cost, node2, cost_pointer_w, cost2[costi(parent2)]
                    )

                cost_pointer_vw.path = strategy_path
                strategy_v[w_in_pre_ltr] = cost_pointer_vw.path

            if node1.children:
                for cost in cost_pointer_v:
                    cost.set_lri(0)
                rows_to_reuse.append(cost_pointer_v)

        return strategy

    def gted(self):
        """Implements GTED algorithm [1, Section 3.4].
        TODO: Document the internals. Point to lines of the algorithm.

        Return the tree edit distance between the source and destination trees.
        """
        it1, it2 = self.it1, self.it2
        subtree1, subtree2 = it1.current_node, it2.current_node

        if subtree1.size == 1 or subtree2.size == 1:  # Use spf1
            return spf1(it1, it2, self.config, subtree1, subtree2)

        path_id = int(self.delta[subtree1.pre_ltr][subtree2.pre_ltr])
        node_id = abs(path_id) - 1

        if node_id < it1.tree_size:  # Apply on subtree 1
            return self.sub_gted(it1, it2, subtree1, path_id, node_id, False)

        # Apply on subtree 2
        node_id -= it1.tree_size
        return self.sub_gted(it2, it1, subtree2, path_id, node_id, True)

    def sub_gted(self, it_f, it_s, subtree_f, path_id, node_id, reverse=False):
        """Apply gted to subtree"""
        # pylint: disable=too-many-arguments
        size = self.it1.tree_size

        strategy = self.get_strategy_path_type(path_id, size, subtree_f)
        current = it_f.pre_ltr_info[node_id]

        for parent, last in it_f.traverse_up(current, subtree_f):
            for child in parent.children:
                if child is not last:
                    it_f.current_index = child.pre_ltr
                    self.gted()

        # todo: Move this property away from node indexer and pass
        # directly to spfs.
        it_f.current_index = subtree_f.pre_ltr

        # Pass to spfs a boolean that says says if the order of input
        # subtrees has been swapped compared to the order of the initial
        # input trees. Used for accessing delta array and deciding on the
        # edit operation. See [1, Section 3.4].

        return SinglePathFunction(
            it_f, it_s, self, node_id, strategy, reverse
        )()

    def get_strategy_path_type(self, strategy_path_id, size1, current):
        """Decodes the path from the optimal strategy to its type.

        Params:
          strategy_path_id: raw path id from strategy array.
          size1: offset used to distinguish between paths in the source and
            destination trees.
          it: node indexer
          current: current subtree processed in tree decomposition phase

        Return type of the strategy path: LEFT, RIGHT, INNER"""
        # pylint: disable=no-self-use
        if math.copysign(1, strategy_path_id) == -1:
            return LEFT
        path_id = abs(strategy_path_id) - 1
        if path_id >= size1:
            path_id = path_id - size1
        if path_id == (current.pre_ltr + current.size - 1):
            return RIGHT
        return INNER

    def compute_edit_mapping(self):
        """Compute the edit mapping between two trees. The trees are input trees
        to the distance computation and the distance must be computed before
        computing the edit mapping (distances of subtree pairs are required)

        Returns list of pairs of nodes that are mapped as pairs
        Nodes that are delete or inserted are mapped to 0
        """
        # todo: Mapping computation requires more thorough documentation
        #       (methods computeEditMapping, forestDist, mappingCost).
        # todo: Before computing the mapping, verify if TED has been computed
        #       Mapping computation should trigger distance computation if
        #       necessary
        self.compute_edit_distance()
        post_ltr_1, post_ltr_2 = self.it1.post_ltr_info, self.it2.post_ltr_info
        size1, size2 = self.it1.tree_size, self.it2.tree_size
        delete, insert = self.config.delete, self.config.insert

        forestdist = zeros((size1 + 1, size2 + 1), float)
        root_node_pair = True
        # Forestdist for input trees has to be computed
        self.forest_dist(size1, size2, forestdist)

        # Empty edit mapping
        edit_mapping = []
        # Stack of tree pairs starting with the pair (ted1, ted2)
        tree_pairs = [(size1, size2)]

        while tree_pairs:
            # Get next pair to be processed
            row, col = tree_pairs.pop()

            # compute forest distance matrix
            if not root_node_pair:
                self.forest_dist(row, col, forestdist)

            root_node_pair = False

            # compute mapping for current forest distance matrix
            first_row = post_ltr_1[row - 1].lld.post_ltr
            first_col = post_ltr_2[col - 1].lld.post_ltr

            while row > first_row or col > first_col:
                row_m1 = post_ltr_1[row - 1]
                col_m1 = post_ltr_2[col - 1]
                fdrc = forestdist[row][col]

                if (
                        row > first_row and
                        forestdist[row - 1][col] + delete(row_m1.node) == fdrc
                ):
                    # Node with post ltr row is deleted from ted1
                    edit_mapping.append((row, 0))
                    row -= 1
                elif (
                        col > first_col and
                        forestdist[row][col - 1] + insert(col_m1.node) == fdrc
                ):
                    # Node with post ltr col is inserted intro ted2
                    edit_mapping.append((0, col))
                    col -= 1
                else:
                    # Node with post ltr row in ted1 is renamed to node col in
                    # ted 2
                    row_lld = row_m1.lld.post_ltr
                    col_lld = col_m1.lld.post_ltr
                    if row_lld == first_row and col_lld == first_col:
                        edit_mapping.append((row, col))
                        row -= 1
                        col -= 1
                    else:
                        # append subtree pair
                        tree_pairs.append((row, col))
                        row, col = row_lld, col_lld

        return edit_mapping


    def forest_dist(self, i, j, forestdist):
        """Recalculates distances between subforests of two subtrees.
        These values are used in mapping computation to track back the origin of
        minimum values. It is basen on Zhang and Shasha algorithm.

        The rename cost must be added in the last line. Otherwise the formula is
        incorrect. This is due to delta storing distances between subtrees
        without the root nodes.

        i and j are postorder ids of the nodes - starting with 1.
        """
        delete, insert = self.config.delete, self.config.insert
        rename, delta = self.config.rename, self.delta

        post_ltr_1, post_ltr_2 = self.it1.post_ltr_info, self.it2.post_ltr_info

        i_m1, j_m1 = post_ltr_1[i - 1], post_ltr_2[j - 1]
        i_m1_lld, j_m1_lld = i_m1.lld.post_ltr, j_m1.lld.post_ltr
        forestdist[i_m1_lld][j_m1_lld] = 0

        for di in range(i_m1_lld + 1, i + 1):
            di_m1 = post_ltr_1[di - 1]
            forestdist[di][j_m1_lld] = (
                forestdist[di - 1][j_m1_lld] + delete(di_m1.node)
            )
            for dj in range(j_m1_lld + 1, j + 1):
                dj_m1 = post_ltr_2[dj - 1]
                forestdist[i_m1_lld][dj] = (
                    forestdist[i_m1_lld][dj - 1] + insert(dj_m1.node)
                )
                cost_ren = rename(di_m1.node, dj_m1.node)
                # todo: the first two elements of the minimum can be computed
                # here, similarly to spfl and spfr
                if di_m1.lld.post_ltr == i_m1_lld and dj_m1.lld.post_ltr == j_m1_lld:
                    forestdist[di][dj] = min(
                        forestdist[di - 1][dj] + delete(di_m1.node),
                        forestdist[di][dj - 1] + insert(dj_m1.node),
                        forestdist[di - 1][dj - 1] + cost_ren
                    )
                    # If substituted with delta, this will overwrite the value
                    # in delta.
                    # It looks that we don't have to write this value.
                    # Conceptually it is correct because we already have all
                    # the values in delta for subtrees without the root nodes,
                    # and we need these.
                    # treedist[di][dj] = forestdist[di][dj];
                else:
                    # di and dj are postorder ids of the nodes - starting with 1
                    # Substituted 'treedist[di][dj]' with
                    # 'delta[it1.postL_to_preL[di-1]][it2.postL_to_preL[dj-1]]'
                    forestdist[di][dj] = min(
                        forestdist[di - 1][dj] + delete(di_m1.node),
                        forestdist[di][dj - 1] + insert(dj_m1.node),
                        forestdist[di_m1.lld.post_ltr][dj_m1.lld.post_ltr] +
                        delta[di_m1.pre_ltr][dj_m1.pre_ltr] + cost_ren
                    )

    def mapping_cost(self, mapping):
        """Calculates the cost of an edit mapping. It traverses the mapping and
        sums up the cost of each operation. The costs are taken from the cost
        model."""
        post_ltr_1, post_ltr_2 = self.it1.post_ltr_info, self.it2.post_ltr_info
        delete, insert = self.config.delete, self.config.insert
        rename = self.config.rename
        cost = 0
        for row, col in mapping:
            if row == 0: # insertion
                cost += insert(post_ltr_2[col - 1].node)
            elif col == 0: # deletion
                cost += delete(post_ltr_1[row - 1].node)
            else:
                cost += rename(post_ltr_1[row - 1].node, post_ltr_2[col - 1].node)
        return cost
