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
"""Helpers for APTED. You can have your own version of theses classes with the
same interface"""

# pylint: disable=no-self-use
# pylint: disable=unused-argument
from __future__ import (absolute_import, division)

from collections import deque
from itertools import chain


class Config(object):
    """Algorithm configuration"""

    valuecls = int

    def delete(self, node):
        """Calculates the cost of deleting a node"""
        return 1

    def insert(self, node):
        """Calculates the cost of inserting a node"""
        return 1

    def rename(self, node1, node2):
        """Calculates the cost of renaming the label of the source node
        to the label of the destination node"""
        return int(node1.name != node2.name)

    def children(self, node):
        """Returns children of node"""
        return getattr(node, 'children', [])

    def mapping_cost(self, apted, mapping):
        """Calculates the cost of an edit mapping. It traverses the mapping and
        sums up the cost of each operation. The costs are taken from the cost
        model."""
        delete, insert = self.delete, self.insert
        rename = self.rename
        cost = 0
        for row, col in mapping:
            if row is None: # insertion
                cost += insert(col.node)
            elif col is None: # deletion
                cost += delete(row.node)
            else:
                cost += rename(row.node, col.node)
        return cost



class PerEditOperationConfig(Config):
    """Algorithm configuration"""

    valuecls = int

    def __init__(self, del_cost, ins_cost, ren_cost):
        self.del_cost = del_cost
        self.ins_cost = ins_cost
        self.ren_cost = ren_cost

    def delete(self, node):
        """Calculates the cost of deleting a node"""
        return self.del_cost

    def insert(self, node):
        """Calculates the cost of inserting a node"""
        return self.ins_cost

    def rename(self, node1, node2):
        """Calculates the cost of renaming the label of the source node
        to the label of the destination node"""
        return (
            super(PerEditOperationConfig, self).rename(node1, node2) *
            self.ren_cost
        )


class ChainedValue(object):
    """Represents a chained value"""
    # pylint: disable=too-few-public-methods

    def __init__(self, value=0, chain_=None):
        self.value = value
        self.chain = chain_ or []

    def __add__(self, other):
        return ChainedValue(
            self.value + other.value, chain(self.chain, other.chain)
        )

    def __sub__(self, other):
        return ChainedValue(
            self.value - other.value, chain(self.chain, [("R", other.chain)])
        )

    def __neg__(self):
        return ChainedValue(
            -self.value, [("R", self.chain)]
        )

    def __bool__(self):
        return bool(self.value)

    def __nonzero__(self):
        return bool(self.value)

    def __repr__(self):
        return repr(self.value)

    def __str__(self):
        return str(self.value)

    def __hash__(self):
        return hash(self.value)

    def __lt__(self, other):
        return self.value < other.value

    def __eq__(self, other):
        return self.value == other.value

    def __int__(self):
        return self.value


def meta_chained_config(config_cls):
    """Creates a config class that keeps track of the chain"""
    class ChainedConfig(config_cls):
        """Chained config class"""

        valuecls = ChainedValue

        def delete(self, node):
            """Calculates the cost of deleting a node"""
            return ChainedValue(
                super(ChainedConfig, self).delete(node),
                [(node, None)]
            )

        def insert(self, node):
            """Calculates the cost of inserting a node"""
            return ChainedValue(
                super(ChainedConfig, self).insert(node),
                [(None, node)]
            )

        def rename(self, node1, node2):
            """Calculates the cost of renaming the label of the source node
            to the label of the destination node"""
            return ChainedValue(
                super(ChainedConfig, self).rename(node1, node2),
                [(node1, node2)]
            )

        def compute_edit_mapping(self, apted):
            """Compute the edit mapping between two trees.

            Returns list of pairs of nodes that are mapped as pairs
            Nodes that are delete or inserted are mapped to None
            """
            value = apted.compute_edit_distance()
            if apted.mapping is None:
                result = set()
                rem_list = set()
                for pair in value.chain:
                    if pair[0] == "R":
                        for rem_pair in pair[1]:
                            try:
                                result.remove(rem_pair)
                            except KeyError:
                                rem_list.add(rem_pair)
                    elif pair not in rem_list:
                        result.add(pair)
                    else:
                        rem_list.remove(pair)
                apted.mapping = result
            return apted.mapping

        def mapping_cost(self, apted, mapping):
            """Calculates the cost of an edit mapping. It traverses the mapping and
            sums up the cost of each operation. The costs are taken from the cost
            model."""
            cost = 0
            for row, col in mapping:
                if row is None: # insertion
                    cost += config_cls.insert(self, col)
                elif col is None: # deletion
                    cost += config_cls.delete(self, row)
                else:
                    cost += config_cls.rename(self, row, col)
            return cost

    return ChainedConfig


class Tree(object):
    """Represents a Tree Node"""

    def __init__(self, name, *children):
        self.name = name
        self.children = list(children)

    def bracket(self):
        """Show tree using brackets notation"""
        result = str(self.name)
        for child in self.children:
            result += child.bracket()
        return "{{{}}}".format(result)

    def __repr__(self):
        return self.bracket()

    @classmethod
    def from_text(cls, text):
        """Create tree from bracket notation

        Bracket notation encodes the trees with nested parentheses, for example,
        in tree {A{B{X}{Y}{F}}{C}} the root node has label A and two children
        with labels B and C. Node with label B has three children with labels
        X, Y, F.
        """
        tree_stack = []
        stack = []
        for letter in text:
            if letter == "{":
                stack.append("")
            elif letter == "}":
                text = stack.pop()
                children = deque()
                while tree_stack and tree_stack[-1][1] > len(stack):
                    child, _ = tree_stack.pop()
                    children.appendleft(child)

                tree_stack.append((cls(text, *children), len(stack)))
            else:
                stack[-1] += letter
        return tree_stack[0][0]
