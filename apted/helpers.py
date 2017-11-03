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

class Config(object):
    """Algorithm configuration"""

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


class PerEditOperationConfig(Config):
    """Algorithm configuration"""

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
