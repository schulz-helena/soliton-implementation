"""Soliton automaton.
"""
import copy

import networkx as nx
import numpy as np
from soliton_classes.soliton_graph import SolitonGraph
from soliton_classes.soliton_path import SolitonPath


class SolitonAutomaton:
    """Representation of a soliton automaton, which contains all soliton paths between all pairs of exterior nodes.
    """

    def __init__(self, soliton_graph: SolitonGraph):
        """Initializes a `SolitonAutomaton` object by using a soliton graph and a stop number.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Soliton graph the automaton is based on."""
        self.deterministic: bool
        """Whether the soliton automaton is deterministic."""
        self.strongly_deterministic: bool
        """Whether the soliton automaton is strongly deterministic."""
        self.reachability_deterministic: bool
        """Whether the soliton automaton is reachability-deterministic."""
        self.degree_of_nondeterminism: int
        """The degree of non-determinism of the soliton automaton"""
        self.states_plus_soliton_paths: dict
        """All states of the soliton automaton plus all soliton paths that can be found in each state
        (Id/ string of the states adjacency matrix as key and state and soliton paths as values)."""
        self.deterministic, self.strongly_deterministic, self.reachability_deterministic, self.degree_of_nondeterminism, self.states_plus_soliton_paths = self.all_paths_and_determinism()


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
            if isinstance(path[0], list): # if path is not a real soliton path
                soliton_path = SolitonPath(soliton_graph, path[0])
                soliton_paths.append([soliton_path, path[1]])
            else:
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

    
    def find_all_paths_given_nodes(self, graph: nx.Graph, bindings: dict, end: int, path: list, akt: int, bind: int, paths: list, bindings_all_timesteps: list, poss_sucs_all_timesteps: list):
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
            bindings_all_timesteps (list): Binding types of all edges in the graph for all past timesteps.
            poss_sucs_all_timesteps (list): Possible successor nodes for all past timesteps.

        Returns:
            list: All found paths (is empty if no path exists).
        """
        possible_suc_nodes = [node for node in list(nx.neighbors(graph, akt)) if node != path[len(path)-2] and bindings[tuple(sorted((akt, node)))] != bind and (node == end or node not in self.soliton_graph.exterior_nodes)] # soliton is not allowed to make a direct turnaround and edge to next node has to have the right binding type
        # base case 1: if end node is reached then finished path to paths
        if end == akt and len(path) >= 3:
            paths.append(path)
            return paths

        # base case 2: if we have two successor-equivalent configurations in the configuration trail
        count = 1
        for k in range(0, len(bindings_all_timesteps)-1):
            if bindings == bindings_all_timesteps[k] and akt == path[k] and possible_suc_nodes == poss_sucs_all_timesteps[k]: # if we already had that exact graph, position and successor positions in this configuration trail
                count += 1
                if count == 2:
                    paths.append([path[:k+1], k]) # save path up to the first of the two successor-equivalent configurations and save index k
                    return paths

        # iterate over all nodes that are adjacent to latest node in path
        for node in possible_suc_nodes:
            akt_copy, path_copy, bindings_copy, bind_copy = self.build_copies(akt, path, bindings, bind) # make copies so we can backtrack later
            path.append(node)
            bind = bindings[tuple(sorted((akt, node)))] # change bind to binding type of edge that was just traversed
            bindings = self.change_bindings(bindings, tuple(sorted((akt, node))))
            bindings_copy2 = copy.deepcopy(bindings)
            akt = node
            bindings_all_timesteps.append(bindings_copy2)
            poss_sucs_all_timesteps.append(possible_suc_nodes)
            # call function recursively: if we can find a path if we go further with the decision we just made (with the node we just added) then paths is changed (new path added)
            paths = self.find_all_paths_given_nodes(graph, bindings, end, path, akt, bind, paths, bindings_all_timesteps, poss_sucs_all_timesteps)
            # try different decisions (nodes) next, so we need variables in the state before the last decision
            akt = akt_copy
            path = path_copy
            bindings = bindings_copy
            bind = bind_copy
            bindings_all_timesteps.pop()
            poss_sucs_all_timesteps.pop()

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
        bindings_all_timesteps = [soliton_graph_copy.bindings] # already add bindings for timestep 0
        poss_sucs_all_timesteps = []
        paths = self.find_all_paths_given_nodes(graph, bindings, end, path, akt, bind, paths, bindings_all_timesteps, poss_sucs_all_timesteps)
        soliton_paths = self.build_soliton_paths(paths, soliton_graph)
        
        return soliton_paths


    def all_paths_and_determinism(self):
        """Calls `call_find_all_paths_given_nodes` for all states of the automaton and all pairs of exterior nodes in order to get all possible soliton paths in all states.
        Checks for determinism with the help of the found states and soliton paths.

        Returns:
            bool: Whether the soliton automaton is deterministic or not.
            bool: Whether the soliton automaton is strongly deterministic or not.
            bool: Whether the soliton automaton is reachability-deterministic or not.
            int: The degree of non-determinism of the soliton automaton.
            dict: All states of the soliton graph plus all soliton paths that can be found in each state.
        """
        ext_nodes = []
        initial_matrix = nx.to_numpy_array(self.soliton_graph.graph)
        states = [self.soliton_graph] # stores all possible states as soliton graphs, needed to iterate over states 
        states_plus_soliton_paths = {self.matrix_to_string(initial_matrix): [self.soliton_graph, []]} # stores all states plus all soliton paths that can be found in that state
        deterministic = True
        strongly_deterministic = True
        reachability_deterministic = True
        degree_of_nondeterminism = 1
        state_sucstate_pair = dict() # dictionary with number of found paths for each state + successor state + pair of exterior nodes
        state_sucstate = dict() # dictionary with bool for each state + successor state (first set to False and then only set to True if there is a burst such that there is only one path between these states)
        for key in self.soliton_graph.exterior_nodes:
            ext_nodes.append(key)
        for state in states: # for all states/ soliton graphs of the automaton
            state_matrix_id = self.matrix_to_string(nx.to_numpy_array(state.graph))
            all_paths = []
            for i in range(0, len(ext_nodes)):
                for j in range(0, len(ext_nodes)): # loop over all pairs of exterior nodes 
                    suc_states_matrix_ids = []
                    paths = self.call_find_all_paths_given_nodes(ext_nodes[i], ext_nodes[j], state) # find soliton paths with all pairs of exterior nodes
                    first_real_path = -1
                    maybe_imperf = []
                    real_paths_this_nodepair = []
                    for p, path in enumerate(paths):
                        if isinstance(path, SolitonPath): # if path is a real path
                            all_paths.append(path)
                            real_paths_this_nodepair.append(path)
                            if first_real_path == -1:
                                first_real_path = p
                            resulting_matrix = path.adjacency_matrices_list[len(path.adjacency_matrices_list)-1] # adjacency matrix of the soliton graph the path results in
                            resulting_matrix_id = self.matrix_to_string(resulting_matrix)
                            if resulting_matrix_id not in states_plus_soliton_paths: # if we found a new state
                                states_plus_soliton_paths[resulting_matrix_id] = [path.resulting_soliton_graph, []] # add it to the dictionary
                                states.append(path.resulting_soliton_graph)
                            if resulting_matrix_id not in suc_states_matrix_ids: # if we found a new state
                                suc_states_matrix_ids.append(resulting_matrix_id)
                                if len(suc_states_matrix_ids) > degree_of_nondeterminism:
                                    degree_of_nondeterminism = len(suc_states_matrix_ids)
                            if np.array_equal(resulting_matrix, paths[first_real_path].adjacency_matrices_list[len(paths[first_real_path].adjacency_matrices_list)-1]) == False:
                                deterministic = False # two or more different soliton graphs have emerged although the same pair of exterior nodes was used
                                reachability_deterministic = False
                            if (state_matrix_id, resulting_matrix_id, (i,j)) not in state_sucstate_pair:
                                state_sucstate_pair[(state_matrix_id, resulting_matrix_id, (i,j))] = 1
                            else: state_sucstate_pair[(state_matrix_id, resulting_matrix_id, (i,j))] += 1
                            state_sucstate[(state_matrix_id, resulting_matrix_id)] = False # put in an entry in the dictionary for these two states
                        else:
                            maybe_imperf.append([path, state_matrix_id, (i,j)]) # this uncompleted path between these two states and these two exterior nodes may be an imperfect path (or no path at all)
                    imperf_found, state_sucstate_pair = self.identify_imperf_paths(real_paths_this_nodepair, maybe_imperf, state_sucstate_pair)
                    if (len(paths) - len(maybe_imperf)) > 1:
                        strongly_deterministic = False # more than one soliton path was found between one pair of exterior nodes
                    elif (len(paths) - len(maybe_imperf)) == 1 and imperf_found:
                        strongly_deterministic = False # we have only one perfect soliton path but at least one imperfect soliton path as well
            states_plus_soliton_paths[state_matrix_id][1] = all_paths # add all found soliton paths in this state
        
        for key in state_sucstate_pair:
            if state_sucstate_pair[key] == 1:
                state_sucstate[key[0], key[1]] = True # there is a pair of exterior nodes such that these is only one path between this state and this successor state
        for key in state_sucstate:
            if state_sucstate[key] == False: # no such pair of exterior nodes exists for these two states
                reachability_deterministic = False
                break


        return deterministic, strongly_deterministic, reachability_deterministic, degree_of_nondeterminism, states_plus_soliton_paths
    

    def identify_imperf_paths(self, real_paths: list, maybe_imperf: list, state_sucstate_pair: dict):
        """Uses the list of perfect paths and a list of possible imperfect paths for a state and a pair of nodes to identify which paths are actually imperfect paths. 

        Args:
            real_paths (list): All real/ perfect paths that were found.
            maybe_imperf (list): All paths that might be imperfect paths.
            state_sucstate_pair (dict): Dictionary that contains a triple of state, successor state and a pair of nodes as key and the number of found paths between the two states with the pair of nodes as value.

        Returns:
            bool: Whether an imperfect path was found.
            dict: Modified `state_sucstate_pair`; if imperfect path was found then number of found paths for state, successor state and pair of nodes was raised.
        """
        imperf_found = False
        for candidate in maybe_imperf:
            for path in real_paths:
                if candidate[0][0].path_labels == path.path_labels[0:len(candidate[0][0].path_labels)]: # [0][0] since only the first element in an element in maybe_imperf is the candidate and only the first element of the candidate is a path 
                    imperf_found = True # we found an imperfect path
                    resulting_matrix = path.adjacency_matrices_list[len(path.adjacency_matrices_list)-1] # adjacency matrix of the soliton graph the path results in
                    resulting_matrix_id = self.matrix_to_string(resulting_matrix)
                    state_sucstate_pair[(candidate[1], resulting_matrix_id, candidate[2])] += 1 # an imperfect path between these two states and with this pair of exterior nodes was found
        return imperf_found, state_sucstate_pair


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
            if isinstance(soliton_path, SolitonPath):
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
    