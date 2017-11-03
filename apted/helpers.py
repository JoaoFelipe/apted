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

    def show(self, ident=0):
        """Show tree using vertical space"""
        result = "." * ident + str(self.name)
        for child in self.children:
            result += "\n" + child.show(ident + 2)
        return result

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
        """Create tree from bracket notation"""
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
