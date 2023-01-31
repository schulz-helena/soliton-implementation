"""Soliton automata.
"""
import copy
import re

import networkx as nx
import numpy as np
from soliton_classes.soliton_graph import SolitonGraph
from soliton_classes.traversal import Traversal


class MultiwaveSolitonAutomata:
    """Representation of a multiwave soliton automata, which finds all traversals for a burst.
    """

    def __init__(self, soliton_graph: SolitonGraph, bursts: str, stop: int):
        """Initializes a `MultiwaveSolitonAutomata` object by using a soliton graph and a burst.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Soliton graph the automata is based on."""
        self.bursts = bursts
        """The set of input burst."""
        self.stop = stop
        """After how many equal soliton graph + soliton positions to stop searching for traversals on current path in search tree."""
        self.bursts_dicts: list = self.build_bursts_dicts()
        """List of the bursts as dictionaries (soliton number as key and a list containing exterior nodes and entry time as value)."""
        self.deterministic: bool
        """Whether the multiwave soliton automata is deterministic."""
        self.strongly_deterministic: bool
        """Whether the multiwave soliton graph is strongly deterministic."""
        self.states_plus_traversals: dict
        """All states of the soliton automata plus all traversals that can be found in each state plus number of traversals found for each burst
        (Id/ string of the states adjacency matrix as key and state, traversals and list of numbers as values)."""
        self.deterministic, self.strongly_deterministic, self.states_plus_traversals = self.all_traversals()


    def build_bursts_dicts(self):
        """Builds a list of burst dictionaries.

        Returns:
            list: The computed list of burst dictionaries.
        """
        bursts = self.bursts.split(";")
        bursts_dicts = []
        for burst in bursts:
            burst = re.sub(r"[{}]+", "", burst)
            burst_dict = self.burst_dict(burst)
            bursts_dicts.append(burst_dict)

        return bursts_dicts

    def burst_dict(self, burst: str):
        """Builds a dictionary as an internal representation of a burst.

        Args:
            burst (str): The burst string that should be turned into a dictionary.

        Returns:
            dict: The computed dictionary with soliton number as key and a list containing exterior nodes and entry time as value.
        """
        burst_copy = copy.copy(burst)
        burst_dict = {}
        entry = 0 # entry time of the current soliton, is increased every time there is another soliton with a k_i
        soliton = 1 # number of the current soliton
        first = re.search(r"[(][0-9]+[,][0-9]+[)]", burst_copy) # first soliton in burst
        nodes = re.findall(r"[0-9]+", first.group()) # find the two exterior nodes
        burst_dict[soliton] = [self.soliton_graph.exterior_nodes_reverse[nodes[0]], self.soliton_graph.exterior_nodes_reverse[nodes[1]], 0] # first soliton entries at time step 0 (take node ids, not node labels)
        burst_copy = burst_copy[first.span()[1]:] # remove all the information on the first soliton from the burst
        while burst_copy != "":
            found = re.search(r"[||][0-9]+[(][0-9]+[,][0-9]+[)]", burst_copy)
            all_nums = re.findall(r"[0-9]+", found.group())
            entry = entry + int(all_nums[0]) # add the current k_i to the entry time
            soliton += 1
            burst_dict[soliton] = [self.soliton_graph.exterior_nodes_reverse[all_nums[1]], self.soliton_graph.exterior_nodes_reverse[all_nums[2]], entry] # take node ids, not node labels
            burst_copy = burst_copy[found.span()[1]:] # remove from burst
            
        return burst_dict


    def build_traversals(self, travs: list, soliton_graph: SolitonGraph):
        """Turns traversals into objects of class `Traversal`.

        Args:
            travs (list): Traversals that are represented as soliton positions and bindings for each timestep.
            soliton_graph (SolitonGraph): Soliton graph the traversal was found in.

        Returns:
            list: Traversals.
        """
        traversals = []
        for trav in travs:
            if isinstance(trav[len(trav)-1], int): # if trav is not a real traversal but part of an endlessly looping traversal
                traversal = Traversal(soliton_graph, trav[0])
                traversals.append([traversal, trav[1]])
            else:
                traversal = Traversal(soliton_graph, trav)
                traversals.append(traversal)

        return traversals


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

    
    def find_all_travs_given_burst (self, graph: nx.Graph, burst_dict: dict, t: int, s_pos_all_timesteps: list, bindings_all_timesteps: list, binds_all_timesteps: list, travs: list):
        """Finds all possible traversals for a given burst by using a recursive backtracking algorithm.

        Args:
            graph (nx.Graph): Graph the traversals should be found in.
            burst_dict (dict): The burst that is used. 
            t (int): Current timestep.
            s_pos_all_timesteps (list): Current soliton positions for each timestep.
            bindings_all_timesteps (list): Current binding types of all edges in the graph.
            binds_all_timesteps (list): Binding type of the last edge that was traversed by each soliton.
            travs (list): All currently found traversals.

        Returns:
            list: All found traversals (is empty if no traversal exists).
        """
        t += 1
        akt_bindings = bindings_all_timesteps[t-1]
        akt_binds = binds_all_timesteps[t-1]
        akt_positions = s_pos_all_timesteps[t-1]

        finished = True
        for soliton in akt_positions:
            if (akt_positions[soliton] != burst_dict[soliton][1] and akt_positions[soliton] != -2) or s_pos_all_timesteps[t-2][soliton] == -1 or t < 2: # if not all solitons reached their end node or already left the graph
                finished = False
                break
        # base case 1: if all solitons traversed the graph successfully
        if finished == True:
            trav = []
            for i, pos in enumerate(s_pos_all_timesteps):
                this_timestep = (pos, bindings_all_timesteps[i]) # for each timestep, add a tuple of the solitons' positions and the bindings
                trav.append(this_timestep)
            travs.append(trav) # trav = traversal of all solitons through the graph
            return travs

        # base case 2: if some solitons are stuck in an endless loop
        count = 1
        for k in range(0, len(bindings_all_timesteps)-1):
            if akt_bindings == bindings_all_timesteps[k] and akt_positions == s_pos_all_timesteps[k]: # if we already had that exact graph and position map in this configuration trail
                count += 1
                if count == self.stop:
                    trav = []
                    for i, pos in enumerate(s_pos_all_timesteps):
                        this_timestep = (pos, bindings_all_timesteps[i])
                        trav.append(this_timestep)
                    travs.append([trav, k+1]) # append the found trav plus the loop point/ timestep (+1, because otherwise in animation we would display the loop point twice)
                    return travs
        

        # in this timestep find possible edges for each soliton individually
        possible_edges = {} # possible edges for each soliton in this timestep
        possible_nodes = {} # possible nodes for each soliton in this timestep
        for soliton in burst_dict:
            end = burst_dict[soliton][1] # end node for this soliton
            akt_pos = akt_positions[soliton] # current position of this soliton
            possible_edges[soliton] = []
            possible_nodes[soliton] = []
            if akt_pos == -2 or (akt_pos == end and s_pos_all_timesteps[t-2][soliton] != -1 and t > 1): # if soliton already traversed the graph successfully
                possible_edges[soliton].append(-2)
                possible_nodes[soliton].append(-2)
            elif burst_dict[soliton][2] > t: # if soliton's entry time has not come yet
                possible_edges[soliton].append(-1)
                possible_nodes[soliton].append(-1)
            elif burst_dict[soliton][2] == t: # if soliton enters the graph now (place it at its starting node)
                possible_edges[soliton].append(burst_dict[soliton][0])
                possible_nodes[soliton].append(burst_dict[soliton][0])
            else:
                for node in list(nx.neighbors(graph, akt_pos)): # iterate over all nodes that are adjacent to current position
                    # soliton can not go to an exterior node that is not the end node, soliton is not allowed to make a direct turnaround and edge to next node has to have the right binding type
                    if ((node not in self.soliton_graph.exterior_nodes) or (node == end)) and node != s_pos_all_timesteps[t-2][soliton] and akt_bindings[tuple(sorted((akt_pos, node)))] != akt_binds[soliton]:
                        possible_edges[soliton].append(tuple(sorted((akt_pos, node))))
                        possible_nodes[soliton].append(node)
        # find all possible combinations of the found edges 
        edges_combs = []
        nodes_combs = []
        for soliton in possible_edges:
            if possible_edges[soliton] == []: # if for one soliton there are no possible next edges, no combinations are possible
                edges_combs.clear()
                nodes_combs.clear()
            if soliton != 1 and edges_combs == []: # if no possible edge combination exists anymore, terminate the loop
                break
            for i, edge in enumerate(possible_edges[soliton]): # for all possible edges for this soliton
                if soliton == 1:
                    edges_combs.append([edge]) # simply add the edge if it's the first soliton 
                    nodes_combs.append({soliton: possible_nodes[soliton][i]}) # add dictionary with node for this soliton
                else:
                    for j, comb in enumerate(edges_combs): # if it's not the first soliton loop over all existing combinations
                        same_ext_node = False
                        for sol in nodes_combs[j]: # check if in this case two solitons would enter/ leave the graph via the same exterior node
                            if nodes_combs[j][sol] == possible_nodes[soliton][i] and nodes_combs[j][sol] in self.soliton_graph.exterior_nodes:
                                same_ext_node = True # this is not allowed to happen so combination is invalid
                        if edge == -1 or edge == -2 or (edge not in comb and same_ext_node == False): # if edge is not already in combination (or soliton is outside the graph)
                            if len(comb) < soliton:
                                edges_combs[j].append(edge) # add the edge if this combination doesn't contain an edge for this soliton already
                                nodes_combs[j][soliton] = possible_nodes[soliton][i]
                            else: # else copy the combination, remove the last edge (the edge for this soliton), add the edge and add this new combination
                                comb_copy = copy.deepcopy(comb)
                                comb_copy.pop()
                                comb_copy.append(edge)
                                edges_combs.append(comb_copy)
                                node_comb_copy = copy.deepcopy(nodes_combs[j])
                                node_comb_copy[soliton] = possible_nodes[soliton][i] # overwrite
                                nodes_combs.append(node_comb_copy)
                        if i == len(possible_edges[soliton])-1: # if we are looking at the last possible edge of this soliton
                            if len(comb) < soliton: # and the combination does not contain an edge for each soliton
                                edges_combs[j] = [] # "delete" the combination
                                nodes_combs[j] = []

        edges_combs = [elem for elem in edges_combs if elem != []]
        nodes_combs = [elem for elem in nodes_combs if elem != []]

        # iterate over all possible combinations
        for i, comb in enumerate(nodes_combs):
            s_pos_all_timesteps.append(comb)
            bindings = copy.deepcopy(akt_bindings)
            binds = copy.deepcopy(akt_binds)
            for j, edge in enumerate(edges_combs[i]): # for all edges in chosen combination
                if isinstance(edge, tuple): # binds and bindings only change if "edge" is actually an edge (if edge is -1, -2 or the start node, no edge is traversed)
                    binds[j+1] = bindings[tuple(sorted(edge))] # change bind to binding type of edge that was just traversed
                    bindings = self.change_bindings(bindings, tuple(sorted(edge))) # change bindings
            bindings_all_timesteps.append(bindings)
            binds_all_timesteps.append(binds)
            # call function recursively: if eventually all solitons reach their end node if we go further with the decision we just made then travs is changed (new traversal added)
            travs = self.find_all_travs_given_burst(graph, burst_dict, t, s_pos_all_timesteps, bindings_all_timesteps, binds_all_timesteps, travs)
            # try different decisions next, so we need variables in the state before the last decision
            s_pos_all_timesteps.pop()
            bindings_all_timesteps.pop()
            binds_all_timesteps.pop()

        return travs # if at some point no other traversal could be added, then all possible traversals are found


    def call_find_all_travs_given_burst (self, burst_dict: dict, soliton_graph: SolitonGraph):
        """Initializes some parameters and then calls `find_all_travs_given_burst` with them.

        Args:
            burst_dict (dict): The burst that is used. 
            soliton_graph (SolitonGraph): Soliton graph the traversals should be found in.

        Returns:
            list: All found travs as traversals (returns empty list when no traversal is found).
        """
        soliton_graph_copy = copy.deepcopy(soliton_graph) # working on copy of graph so no unwanted changes are made
        graph = soliton_graph_copy.graph
        t = 0
        bindings_all_timesteps = [soliton_graph_copy.bindings] # already add bindings for timestep 0
        soliton_positions = {}
        binds = {}
        for soliton in burst_dict:
            binds[soliton] = 0
            if burst_dict[soliton][2] == 0: # if soliton enters the graph at timestep 0
                soliton_positions[soliton] = burst_dict[soliton][0] # starts at its start node
            else:
                soliton_positions[soliton] = -1 # soliton is not in graph at timestep 0
        # already add positions and binds for timestep 0
        s_pos_all_timesteps = [soliton_positions]
        binds_all_timesteps = [binds]
        travs = []

        travs = self.find_all_travs_given_burst(graph, burst_dict, t, s_pos_all_timesteps, bindings_all_timesteps, binds_all_timesteps, travs)
        traversals = self.build_traversals(travs, soliton_graph)

        return traversals


    def all_traversals (self):
        """Calls `call_find_all_travs_given_burst` for all states of the automata and all bursts in order to get all possible traversals in all states.

        Returns:
            dict: All states of the soliton graph plus all traversals that can be found in each state plus number of traversals found for each burst.
        """
        initial_matrix = nx.to_numpy_array(self.soliton_graph.graph)
        states = [self.soliton_graph] # stores all possible states as soliton graphs, needed to iterate over states 
        states_plus_traversals = {self.matrix_to_string(initial_matrix): [self.soliton_graph, [], []]} # stores all states plus all traversals that can be found in that state plus number of traversals found for each burst
        deterministic = True
        strongly_deterministic = True

        for state in states: # for all states/ soliton graphs of the automata
            state_matrix_id = self.matrix_to_string(nx.to_numpy_array(state.graph))
            all_traversals = []
            num_traversals_per_burst = []
            for burst_dict in self.bursts_dicts: # loop over all bursts
                traversals = self.call_find_all_travs_given_burst(burst_dict, state) # find soliton paths with all bursts
                num_traversals_per_burst.append(len(traversals))
                loops_num = 0
                first_real_trav = -1
                for t, traversal in enumerate(traversals):
                    all_traversals.append(traversal)
                    if isinstance(traversal, Traversal): # if traversal is a real traversal and no endless loop
                        if first_real_trav == -1:
                            first_real_trav = t
                        resulting_matrix = traversal.adjacency_matrices_list[len(traversal.adjacency_matrices_list)-1] # adjacency matrix of the soliton graph the traversal results in
                        resulting_matrix_id = self.matrix_to_string(resulting_matrix)
                        if resulting_matrix_id not in states_plus_traversals: # if we found a new state
                            states_plus_traversals[resulting_matrix_id] = [traversal.resulting_soliton_graph, [], []] # add it to the dictionary
                            states.append(traversal.resulting_soliton_graph)
                        if np.array_equal(resulting_matrix, traversals[first_real_trav].adjacency_matrices_list[len(traversals[first_real_trav].adjacency_matrices_list)-1]) == False:
                                deterministic = False # two or more different soliton graphs have emerged although the same pair of exterior nodes was used
                    else: loops_num += 1
                if (len(traversals) - loops_num) > 1:
                        strongly_deterministic = False # more than one soliton path was found between one pair of exterior nodes
            states_plus_traversals[state_matrix_id][1] = all_traversals # add all found traversals in this state
            states_plus_traversals[state_matrix_id][2] = num_traversals_per_burst # add number of traversals found for each burst (in this state)

        return deterministic, strongly_deterministic, states_plus_traversals

        
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
