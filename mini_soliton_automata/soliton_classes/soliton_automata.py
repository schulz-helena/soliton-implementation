"""Soliton automata.
"""
import copy

import networkx as nx
import numpy as np
from soliton_classes.soliton_graph import SolitonGraph
from soliton_classes.soliton_path import SolitonPath


class SolitonAutomata:
    """Representation of a soliton automata, which finds all soliton paths between all pairs of exterior nodes.
    """

    def __init__(self, soliton_graph):
        """Initializes a `SolitonAutomata` object by using a soliton graph.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Soliton graph the automata is based on."""
        self.deterministic: bool
        """Whether the soliton graph is deterministic."""
        self.strongly_deterministic: bool
        """Whether the soliton graph is strongly deterministic."""
        self.states_plus_soliton_paths: dict
        """All states of the soliton graph plus all soliton paths that can be found in each state
        (Id/ string of the states adjacency matrix as key and state and soliton paths as values)."""
        self.deterministic, self.strongly_deterministic, self.states_plus_soliton_paths = self.all_paths_and_determinism()
        initial_matrix_id = self.matrix_to_string(nx.to_numpy_array(self.soliton_graph.graph))
        self.soliton_paths: list = self.states_plus_soliton_paths[initial_matrix_id][1]
        """All found soliton paths in initial soliton graph"""


    def build_soliton_paths(self, paths: list, soliton_graph: SolitonGraph):
        """Turns paths into objects of class `SolitonPath`.

        Args:
            paths (list): Paths that are represented as lists of node ids.
            soliton_graph (SolitonGraph): Soliton graph the paths were found in.

        Returns:
            list: Soliton paths.
        """
        soliton_paths = []
        for path in paths:
            soliton_path = SolitonPath(soliton_graph, path)
            soliton_paths.append(soliton_path)

        return soliton_paths


    def change_bindings(self, bindings: dict, edge: tuple):
        """Changes binding type of an edge (1 -> 2 and 2 -> 1).

        Args:
            bindings (dict): Current binding types of all edges in the graph.
            edge (tuple): Edge whose binding type should be changed.

        Returns:
            dict: Updated binding dictionary.
        """
        if bindings[(edge[0], edge[1])] == 2:
            bindings[(edge[0], edge[1])] = 1
        else:
            bindings[(edge[0], edge[1])] = 2
            
        return bindings


    def build_copies(self, akt: int, path: list, bindings: dict, bind: int):
        """Copies all variables that are changed during `find_all_paths`.

        Args:
            akt (int): Node that was currently added to path.
            path (list): Current found path.
            bindings (dict): Current binding types of all edges in the graph.
            bind (int): Binding type of the last edge that was traversed.

        Returns:
            int: Copy of `akt`.
            list: Copy of `path`.
            dict: Copy of `bindings`.
            int: Copy of `bind`.
        """
        akt_copy = copy.deepcopy(akt)
        path_copy = copy.deepcopy(path)
        bindings_copy = copy.deepcopy(bindings)
        bind_copy = copy.deepcopy(bind)

        return akt_copy, path_copy, bindings_copy, bind_copy


    def find_all_paths_given_nodes(self, graph: nx.Graph, bindings: dict, end: int, path: list, akt: int, bind: int, paths: list):
        """Finds all possible soliton paths between two given exterior nodes by using a recursive backtracking algorithm.
            A path can only be a soliton path if the edges traversed by the soliton have alternating binding types (1,2,1,2,...).

        Args:
            graph (nx.Graph): Graph the paths should be found in.
            bindings (dict): Current binding types of all edges in the graph.
            end (int): End node of path.
            path (list): Current found path.
            akt (int): Node that was currently added to path.
            bind (int): Binding type of the last edge that was traversed.
            paths (list): All currently found paths.

        Returns:
            list: All found paths (is empty if no path exists).
        """
        # base case: if end node is reachable then add end node to path and add finished path to paths
        if end in list(nx.neighbors(graph, akt)) and bindings[tuple(sorted((akt, end)))] != bind and end != path[len(path)-2]: # also in this case a direct turnaround is not allowed
            path.append(end)
            bindings = self.change_bindings(bindings, (akt, end)) # change binding of traversed edge
            paths.append(path)
            return paths

        # iterate over all nodes that are adjacent to latest node in path
        for node in list(nx.neighbors(graph, akt)):
            if node != path[len(path)-2] and bindings[tuple(sorted((akt, node)))] != bind: # soliton is not allowed to make a direct turnaround and edge to next node has to have the right binding type
                akt_copy, path_copy, bindings_copy, bind_copy = self.build_copies(akt, path, bindings, bind) # make copies so we can backtrack later
                path.append(node)
                bind = bindings[tuple(sorted((akt, node)))] # change bind to binding type of edge that was just traversed
                bindings = self.change_bindings(bindings, tuple(sorted((akt, node))))
                akt = node
                # call function recursively: if we can find a path if we go further with the decision we just made (with the node we just added) then paths is changed (new path added)
                paths = self.find_all_paths_given_nodes(graph, bindings, end, path, akt, bind, paths)
                # try different decisions (nodes) next, so we need variables in the state before the last decision
                akt = akt_copy
                path = path_copy
                bindings = bindings_copy
                bind = bind_copy

        return paths # if at some point no new node could be added, then all possible paths are found
    

    def call_find_all_paths_given_nodes(self, start: int, end: int, soliton_graph: SolitonGraph):
        """Initializes some parameters and then calls `find_all_paths_given_nodes` with them.

        Args:
            start (int): Start node of path.
            end (int): End node of path.
            soliton_graph (SolitonGraph): Soliton graph the paths should be found in.

        Returns:
            list: All found paths as soliton paths (returns empty list when no path is found).
        """
        paths = []
        soliton_graph_copy = copy.deepcopy(soliton_graph) # working on copy of graph so no unwanted changes are made
        graph = soliton_graph_copy.graph
        bindings = soliton_graph_copy.bindings
        path = [start]
        akt = start
        bind = 0
        paths = self.find_all_paths_given_nodes(graph, bindings, end, path, akt, bind, paths)
        soliton_paths = self.build_soliton_paths(paths, soliton_graph)
        
        return soliton_paths


    def all_paths_and_determinism(self):
        """Calls `call_find_all_paths_given_nodes` for all states of the automata and all pairs of exterior nodes in order to get all possible soliton paths in all states.
        Checks for determinism with the help of the found states and soliton paths.

        Returns:
            bool: Whether the soliton graph is deterministic or not.
            bool: Whether the soliton graph is strongly deterministic or not.
            dict: All states of the soliton graph plus all soliton paths that can be found in each state.
        """
        ext_nodes = []
        initial_matrix = nx.to_numpy_array(self.soliton_graph.graph)
        states = [self.soliton_graph] # stores all possible states as soliton graphs, needed to iterate over states 
        states_plus_soliton_paths = {self.matrix_to_string(initial_matrix): [self.soliton_graph, []]} # stores all states plus all soliton paths that can be found in that state
        deterministic = True
        strongly_deterministic = True
        for key in self.soliton_graph.exterior_nodes:
            ext_nodes.append(key)
        for state in states: # for all states/ soliton graphs of the automata
            state_matrix_id = self.matrix_to_string(nx.to_numpy_array(state.graph))
            all_paths = []
            for i in range(0, len(ext_nodes)):
                for j in range(0, len(ext_nodes)): # loop over all pairs of exterior nodes 
                    paths = self.call_find_all_paths_given_nodes(ext_nodes[i], ext_nodes[j], state) # find soliton paths with all pairs of exterior nodes
                    for path in paths:
                        all_paths.append(path)
                        resulting_matrix = path.adjacency_matrices_list[len(path.adjacency_matrices_list)-1] # adjacency matrix of the soliton graph the path results in
                        resulting_matrix_id = self.matrix_to_string(resulting_matrix)
                        if resulting_matrix_id not in states_plus_soliton_paths: # if we found a new state
                            states_plus_soliton_paths[resulting_matrix_id] = [path.resulting_soliton_graph, []] # add it to the dictionary
                            states.append(path.resulting_soliton_graph)
                        if np.array_equal(resulting_matrix, paths[0].adjacency_matrices_list[len(paths[0].adjacency_matrices_list)-1]) == False:
                            deterministic = False # two or more different soliton graphs have emerged although the same pair of exterior nodes was used
                    if len(paths) > 1:
                        strongly_deterministic = False # more than one soliton path was found between one pair of exterior nodes
            states_plus_soliton_paths[state_matrix_id][1] = all_paths # add all found soliton paths in this state

        return deterministic, strongly_deterministic, states_plus_soliton_paths

    def matrix_to_string(self, matrix: np.ndarray):
        """Computes an ID for a matrix.

        Args:
            matrix (np.ndarray): Matrix whos ID should be computed.

        Returns:
            str: The computed ID which is a concatination of the matrix's elements in row-wise order.
        """
        matrix_id = ""
        for row in matrix:
            row_id = ''.join(str(int(elem)) for elem in row)
            matrix_id = matrix_id + row_id
        return matrix_id

    def find_impervious_paths(self):
        """Finds all impervious paths of the initial soliton graph.

        Returns:
            list: The found impervious paths as readable user outputs.
        """
        unused_edges = list(self.soliton_graph.graph.edges)
        initial_matrix_id = self.matrix_to_string(nx.to_numpy_array(self.soliton_graph.graph))
        initial_soliton_paths = self.states_plus_soliton_paths[initial_matrix_id][1]
        # remove all edges that are traversed in any soliton path:
        for soliton_path in initial_soliton_paths:
            for i, node in enumerate(soliton_path.path):
                if i < (len(soliton_path.path)-1):
                    if tuple(sorted((node, soliton_path.path[i+1]))) in unused_edges:
                        unused_edges.remove(tuple(sorted((node, soliton_path.path[i+1]))))

        akt_path = []
        imperv_paths = []
        for i, edge in enumerate(unused_edges):
            if edge[0] not in akt_path: akt_path.append(edge[0]) # don't add nodes twice
            if edge[1] not in akt_path: akt_path.append(edge[1])
            if i == (len(unused_edges)-1): # edge is last element
                soli_path = SolitonPath(self.soliton_graph, akt_path)
                imperv_paths.append(soli_path.path_for_user)
            elif edge[1] != unused_edges[i+1][0]: # next edge is not connected to current path
                soli_path = SolitonPath(self.soliton_graph, akt_path)
                imperv_paths.append(soli_path.path_for_user)
                akt_path = [] # from next edge on it is a different path

        return imperv_paths
