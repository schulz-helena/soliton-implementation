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

    def __init__(self, soliton_graph: SolitonGraph, burst: str):
        """Initializes a `MultiwaveSolitonAutomata` object by using a soliton graph and a burst.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Soliton graph the automata is based on."""
        self.burst = burst
        """The input burst."""
        self.burst_dict = self.build_burst_dict()
        """The burst as a dictionary (soliton number as key and a list containing exterior nodes and entry time as value)."""
        self.traversals: list
        """All found traversals"""


    def build_burst_dict(self):
        """Builds a dictionary as an internal representation of a burst.

        Returns:
            dict: The computed dictionary with soliton number as key and a list containing exterior nodes and entry time as value.
        """
        burst_copy = copy.copy(self.burst)
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


    def find_all_travs_given_burst (self, t: int, s_pos_all_timesteps: list, bindings_all_timesteps: list, binds_all_timesteps: list, travs: list):
        """Finds all possible traversals for a given burst by using a recursive backtracking algorithm.

        Args:
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
            if akt_positions[soliton] != self.burst_dict[soliton][1] and akt_positions[soliton] != -2: # if not all solitons reached their end node or already left the graph
                finished = False
                break
        # base case: if all solitons traversed the graph successfully
        if finished == True:
            trav = []
            for i, pos in enumerate(s_pos_all_timesteps):
                this_timestep = (pos, bindings_all_timesteps[i]) # for each timestep, add a tuple of the solitons' positions and the bindings
                trav.append(this_timestep)
            travs.append(trav) # trav = traversal of all solitons through the graph
            return travs
        

        # in this timestep find possible edges for each soliton individually
        possible_edges = {} # possible edges for each soliton in this timestep
        possible_nodes = {} # possible nodes for each soliton in this timestep
        for soliton in self.burst_dict:
            end = self.burst_dict[soliton][1] # end node for this soliton
            akt_pos = akt_positions[soliton] # current position of this soliton
            possible_edges[soliton] = []
            possible_nodes[soliton] = []
            if akt_pos == -2 or akt_pos == end: # if soliton already traversed the graph successfully
                possible_edges[soliton].append(-2)
                possible_nodes[soliton].append(-2)
            elif self.burst_dict[soliton][2] > t: # if soliton's entry time has not come yet
                possible_edges[soliton].append(-1)
                possible_nodes[soliton].append(-1)
            elif self.burst_dict[soliton][2] == t: # if soliton enters the graph now (place it at its starting node)
                possible_edges[soliton].append(self.burst_dict[soliton][0])
                possible_nodes[soliton].append(self.burst_dict[soliton][0])
            elif end in list(nx.neighbors(self.soliton_graph.graph, akt_pos)) and akt_bindings[tuple(sorted((akt_pos, end)))] != akt_binds[soliton] and end != s_pos_all_timesteps[t-2][soliton]: # if end node is reachable
                possible_edges[soliton].append(tuple(sorted((akt_pos, end))))
                possible_nodes[soliton].append(end)
            else:
                for node in list(nx.neighbors(self.soliton_graph.graph, akt_pos)): # iterate over all nodes that are adjacent to current position
                    # soliton is not allowed to make a direct turnaround and edge to next node has to have the right binding type
                    if node not in self.soliton_graph.exterior_nodes and node != s_pos_all_timesteps[t-2][soliton] and akt_bindings[tuple(sorted((akt_pos, node)))] != akt_binds[soliton]:
                        possible_edges[soliton].append(tuple(sorted((akt_pos, node))))
                        possible_nodes[soliton].append(node)
        # find all possible combinations of the found edges 
        edges_combs = []
        nodes_combs = []
        for soliton in possible_edges:
            if soliton != 1 and edges_combs == []: # if no possible edge combination exists anymore, terminate the loop
                break
            for i, edge in enumerate(possible_edges[soliton]): # for all possible edges for this soliton
                if soliton == 1:
                    edges_combs.append([edge]) # simply add the edge if it's the first soliton 
                    nodes_combs.append({soliton: possible_nodes[soliton][i]}) # add dictionary with node for this soliton
                else:
                    for j, comb in enumerate(edges_combs): # if it's not the first soliton loop over all existing combinations
                        if edge == -1 or edge == -2 or edge not in comb: # if edge is not already in combination (or soliton is outside the graph)
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
                                del edges_combs[j] # delete the combination
                                del nodes_combs[j]

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
            travs = self.find_all_travs_given_burst(t, s_pos_all_timesteps, bindings_all_timesteps, binds_all_timesteps, travs)
            # try different decisions next, so we need variables in the state before the last decision
            s_pos_all_timesteps.pop()
            bindings_all_timesteps.pop()
            binds_all_timesteps.pop()

        return travs # if at some point no other traversal could be added, then all possible traversals are found


    def call_find_all_travs_given_burst (self, soliton_graph: SolitonGraph):
        """Initializes some parameters and then calls `find_all_travs_given_burst` with them.

        Args:
            soliton_graph (SolitonGraph): Soliton graph the traversals should be found in.

        Returns:
            list: All found travs as traversals (returns empty list when no traversal is found).
        """
        t = 0
        bindings_all_timesteps = [self.soliton_graph.bindings] # already add bindings for timestep 0
        soliton_positions = {}
        binds = {}
        for soliton in self.burst_dict:
            binds[soliton] = 0
            if self.burst_dict[soliton][2] == 0: # if soliton enters the graph at timestep 0
                soliton_positions[soliton] = self.burst_dict[soliton][0] # starts at its start node
            else:
                soliton_positions[soliton] = -1 # soliton is not in graph at timestep 0
        # already add positions and binds for timestep 0
        s_pos_all_timesteps = [soliton_positions]
        binds_all_timesteps = [binds]
        travs = []

        travs = self.find_all_travs_given_burst(t, s_pos_all_timesteps, bindings_all_timesteps, binds_all_timesteps, travs)
        traversals = self.build_traversals(travs, soliton_graph)

        return traversals

        
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
