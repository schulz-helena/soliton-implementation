import copy

import networkx as nx
import numpy as np

from soliton_graph import SolitonGraph


class SolitonPath:
    """Representation of a soliton path.
    """

    def __init__(self, soliton_graph: SolitonGraph, path: list):
        """Initializes a `SolitonPath` object.
        """
        self.path: list = path
        """Path consisting of node ids."""
        self.soliton_graph: SolitonGraph = soliton_graph
        """According soliton graph."""
        self.path_labels: list
        """Path consisiting of node labels."""
        self.path_for_user: str
        """Representation of the path as a readable user output."""
        self.path_labels, self.path_for_user = self.path_representations()
        self.bindings_list: list = self.find_bindings_for_each_timestep()
        """Bindings for each timestep."""
        self.adjacency_matrices_list: list = self.find_adjacency_matrices_for_each_timestep()
        """Adjacency matrices for each timestep."""

    
    def path_representations(self):
        """Creates two representations of a path: One with node labels and one as a string as a readable user output.

        Returns:
            dict: Path with node labels instead of node ids.
            dict: Representation of the path the user gets as an output.
        """
        path_labels = []
        path_for_user = ""
        for i, node in enumerate(self.path):
            path_labels.append(nx.get_node_attributes(self.soliton_graph.graph, 'label')[node])
            if i == len(self.path) - 1:
                path_for_user = path_for_user + f"{path_labels[i]}"
            else:
                path_for_user = path_for_user + f"{path_labels[i]} - "

        return path_labels, path_for_user


    def find_bindings_for_each_timestep(self):
        """Finds the binding dictionary for each timestep of traversal of a path.

        Returns:
            list: Contains bindings for each timestep.
        """
        soliton_graph_copy = copy.deepcopy(self.soliton_graph) # working on a copy of the soliton graph so no unwanted changes occur
        bindings_list = []
        bindings_copy = soliton_graph_copy.bindings.copy()
        bindings_list.append(bindings_copy)
        # for each node in soliton path: change binding of edge that soliton just traversed:
        for i in range(1, len(self.path)): # starting at 1 because self.path[0] has no predecessor node
            if soliton_graph_copy.bindings[tuple(sorted((self.path[i-1], self.path[i])))] == 2:
                soliton_graph_copy.bindings[tuple(sorted((self.path[i-1], self.path[i])))] = 1
            else:
                soliton_graph_copy.bindings[tuple(sorted((self.path[i-1], self.path[i])))] = 2

            bindings_copy = soliton_graph_copy.bindings.copy()
            bindings_list.append(bindings_copy)

        return bindings_list


    def find_adjacency_matrices_for_each_timestep(self):
        """Finds the adjacency matrix for each timestep of traversal of a path.

        Returns:
            list: Contains `Numpy` matrix for each timestep.
        """
        np.set_printoptions(edgeitems=1000, linewidth=100000) # make sure even large matrices are displayed without new lines inside the matrix
        soliton_graph_copy = copy.deepcopy(self.soliton_graph)
        adjacency_matrices_list = []
        for bindings in self.bindings_list:
            soliton_graph_copy.set_bindings(bindings) # change bindings in graph which changes edge weights which changes adjacency matrix
            matrix = nx.to_numpy_matrix(soliton_graph_copy.graph)
            adjacency_matrices_list.append(matrix)

        return adjacency_matrices_list
    
