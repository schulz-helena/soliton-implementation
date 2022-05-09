"""Class that represents a soliton path

    Attributes: path (consisting of node ids), the according soliton graph, the bindings and adjacency matrices for each timestep
"""
import copy
from dataclasses import dataclass

import networkx as nx

from soliton_graph import SolitonGraph


@dataclass
class SolitonPath:

    path: list
    soliton_graph: SolitonGraph

    def __post_init__(self):
        self.bindings_list = self.find_bindings_for_each_timestep()
        self.adjacency_matrices_list = self.find_adjacency_matrices_for_each_timestep()


    def find_bindings_for_each_timestep(self):
        soliton_graph_copy = copy.deepcopy(self.soliton_graph) #working on a copy of the soliton graph so no unwanted changes occur
        bindings_list = []
        bindings_copy = soliton_graph_copy.bindings.copy()
        bindings_list.append(bindings_copy)
        # for each node in soliton path: change binding of edge that soliton just traversed:
        for i in range(1, len(self.path)): # starting at 1 because self.path[0] has no predecessor node
            if (self.path[i-1], self.path[i]) in soliton_graph_copy.bindings:
                if soliton_graph_copy.bindings[(self.path[i-1], self.path[i])] == 2:
                    soliton_graph_copy.bindings[(self.path[i-1], self.path[i])] = 1
                else:
                    soliton_graph_copy.bindings[(self.path[i-1], self.path[i])] = 2
            else:
                if soliton_graph_copy.bindings[(self.path[i], self.path[i-1])] == 2:
                    soliton_graph_copy.bindings[(self.path[i], self.path[i-1])] = 1
                else:
                    soliton_graph_copy.bindings[(self.path[i], self.path[i-1])] = 2
            bindings_copy = soliton_graph_copy.bindings.copy()
            bindings_list.append(bindings_copy)
        return bindings_list


    def find_adjacency_matrices_for_each_timestep(self):
        soliton_graph_copy = copy.deepcopy(self.soliton_graph)
        adjacency_matrices_list = []
        for i in range(0, len(self.bindings_list)):
            bindings = self.bindings_list[i]
            soliton_graph_copy.set_bindings(bindings) #change bindings in graph which changes edge weights which changes adjacency matrix
            matrix = nx.to_numpy_matrix(soliton_graph_copy.graph)
            adjacency_matrices_list.append(matrix)
        return adjacency_matrices_list


if __name__ == "__main__":
    my_graph = SolitonGraph('CC=CC')
    my_path = SolitonPath([0,1,2,3], my_graph)
    print(my_path.bindings_list)
    print(my_path.adjacency_matrices_list)
