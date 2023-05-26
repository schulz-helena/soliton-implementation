"""Traversal.
"""
import copy
import re

import networkx as nx
import numpy as np
from soliton_classes.soliton_graph import SolitonGraph


class Traversal:
    """Representation of a traversal. Traversal means the paths of all solitons in a burst (alternative representation to configurationtrail).
    """

    def __init__(self, soliton_graph: SolitonGraph, pos_and_bindings: list):
        """Initializes a `Traversal` object.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Soliton graph the traversal was found in."""
        self.pos_and_bindings: list = pos_and_bindings
        """Positions and bindings for each timestep."""
        self.pos: list
        """Soliton positions of all solitons for each timestep."""
        self.bindings_list: list
        """Bindings of the graphs edges for each timestep."""
        self.pos, self.bindings_list = self.split_pos_and_bindings()
        self.soliton_num: int = len(self.pos[0])
        """Number of solitons in this traversal."""
        self.traversal_for_user: list
        """Representation of the traversal as a readable user output."""
        self.traversal_node_ids: list
        """Representation of the traversal as seperate paths as lists of node ids"""
        self.traversal_for_user, self.traversal_node_ids = self.traversal_representation()
        self.adjacency_matrices_list: list
        """Adjacency matrices for each timestep."""
        self.resulting_soliton_graph: SolitonGraph
        """Soliton graph the traversal results in."""
        self.adjacency_matrices_list, self.resulting_soliton_graph = self.adjacency_matrices_and_resulting_soliton_graph()


    def split_pos_and_bindings(self):
        """Splits the list containing the positions and the bindings for each timestep into two seperate lists.

        Returns:
            list: Soliton positions of all solitons for each timestep.
            list: Bindings of the graphs edges for each timestep.
        """
        pos = []
        bindings_list = []
        for p_and_b in self.pos_and_bindings:
            pos.append(p_and_b[0])
            bindings_list.append(p_and_b[1])

        return pos, bindings_list


    def traversal_representation(self):
        """Creates a representation of a traversal that is readable for the user. 
        Paths are marked with the corresponding soliton.
        Uses '-' if the soliton is not currently in the graph at a specific timestep and otherwise 
        node labels to indicate where the soliton is.
        Also creates a simple path with just the node ids for each soliton in the traversal.

        Returns:
            list: Representations of soliton paths for each soliton.
            list: Path with node ids for each soliton.
        """
        representations, paths = [], []
        for i in range (1, self.soliton_num+1):
            representations.append(f"S{i}: ") # mark path with number of the soliton
            paths.append([])
        for j, positions in enumerate(self.pos):
            for soliton in positions:
                if positions[soliton] == -2 or positions[soliton] == -1:
                    representations[soliton-1] += "." # currently not in graph
                else:
                    representations[soliton-1] += self.soliton_graph.labels[positions[soliton]] # add node label
                    paths[soliton-1].append(positions[soliton]) # add node id
                if j != len(self.pos)-1:
                    representations[soliton-1] += " - "

        return representations, paths


    def adjacency_matrices_and_resulting_soliton_graph (self):
        """Finds the adjacency matrix for each timestep of a traversal.

        Returns:
            list: Contains `Numpy` matrix for each timestep.
            SolitonGraph: Soliton graph the traversal results in.
        """
        np.set_printoptions(edgeitems=1000, linewidth=100000) # make sure even large matrices are displayed without new lines inside the matrix
        soliton_graph_copy = copy.deepcopy(self.soliton_graph)
        adjacency_matrices_list = []
        for bindings in self.bindings_list:
            soliton_graph_copy.set_bindings(bindings) # change bindings in graph which changes edge weights which changes adjacency matrix
            matrix = nx.to_numpy_array(soliton_graph_copy.graph)
            adjacency_matrices_list.append(matrix)

        return adjacency_matrices_list, soliton_graph_copy
    

    def ring_passages(self):
        paths = self.traversal_node_ids
        rings_num = len(self.soliton_graph.rings)
        ring_passages = dict()
        for i in range(rings_num): # for all rings
            ring_passages[i] = []
            ring = self.soliton_graph.rings[i]
            ring_backw = ring[::-1]
            for j, path in enumerate(paths): # for all paths in traversal (one path for each soliton)
                ring_passages[i].append([]) # a list of ring passages for this ring and this soliton
                pattern = [str(ring[1:])[1:-1], str(ring_backw[1:])[1:-1]]
                regex = re.compile('|'.join(pattern)) # search for "forward" and "backward" ring (doesn't matter if soliton is going left or right when entering ring)
                res = re.finditer(regex, str(path))
                last_end_index = next(res).span()[1] # end index of first match
                count = 1
                for match in res:
                    if match.span()[0] - last_end_index != 2: # if match is not directly behind last match (the two ring passages did not happen back to back)
                        ring_passages[i][j].append(count) # append count of last found back to back ring passages, start again for current ring passage
                        count = 1
                    else:
                        count += 1
                    last_end_index = match.span()[1]
                ring_passages[i][j].append(count)
        return ring_passages
