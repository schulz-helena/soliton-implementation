import copy

import networkx as nx

from soliton_classes.soliton_graph import SolitonGraph
from soliton_classes.soliton_path import SolitonPath


class MiniSolitonAutomata:
    """Representation of part of a soliton automata, which finds all soliton paths for a given pair of exterior nodes.
    """

    def __init__(self, soliton_graph, start, end):
        """Initializes a `MiniSolitonAutomata` object by using a soliton graph and a pair of exterior nodes, represented as node labels.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Soliton graph the automata is based on."""
        self.start: int = self.soliton_graph.exterior_nodes_reverse[str(start)]
        """First exterior node/ start node of soliton paths."""
        self.end: int = self.soliton_graph.exterior_nodes_reverse[str(end)]
        """Second exterior node/ end node of soliton paths."""
        paths = self.call_find_all_paths()
        self.soliton_paths: list = self.build_soliton_paths(paths)
        """All found soliton paths between given pair of exterior nodes."""


    def build_soliton_paths(self, paths: list):
        """Turns paths into objects of class `SolitonPath`.

        Args:
            paths (list): Paths that are represented as lists of node ids.

        Returns:
            list: Soliton paths.
        """
        soliton_paths = []
        for path in paths:
            soliton_path = SolitonPath(self.soliton_graph, path)
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


    def find_all_paths(self, graph: nx.Graph, bindings: dict, end: int, path: list, akt: int, bind: int, paths: list):
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
                paths = self.find_all_paths(graph, bindings, end, path, akt, bind, paths)
                # try different decisions (nodes) next, so we need variables in the state before the last decision
                akt = akt_copy
                path = path_copy
                bindings = bindings_copy
                bind = bind_copy

        return paths # if at some point no new node could be added, then all possible paths are found
    

    def call_find_all_paths(self):
        """Initializes some parameters and then calls `find_all_paths` with them.

        Returns:
            list: All found paths (returns empty list when no path is found).
        """
        paths = []
        soliton_graph_copy = copy.deepcopy(self.soliton_graph) # working on copy of graph so no unwanted changes are made
        graph = soliton_graph_copy.graph
        bindings = soliton_graph_copy.bindings
        path = [self.start]
        akt = self.start
        bind = 0
        res = self.find_all_paths(graph, bindings, self.end, path, akt, bind, paths)
        
        return res
