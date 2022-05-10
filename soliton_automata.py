"""Class that represents a soliton automata and therefore finds all soliton paths for a given pair of exterior nodes 

    Attributes: soliton graph the automata is based on, start and end node for paths, list of soliton paths (right now only the first found soliton path)
"""
import copy
from dataclasses import dataclass

import networkx as nx

from soliton_graph import SolitonGraph
from soliton_path import SolitonPath


@dataclass
class SolitonAutomata:

    soliton_graph: SolitonGraph
    start: int
    end: int

    def __post_init__(self):
        path, bindings_list = self.call_find_first_path()
        self.soliton_path = SolitonPath(path, self.soliton_graph, bindings_list)


    def change_bindings(self, bindings: dict, edge: tuple):
        """Helping function to change binding type of an edge (1 -> 2 and 2 -> 1)

        Args:
            bindings (dict): current binding types of all edges in the graph
            edge (tuple): edge whose binding type should be changed

        Returns:
            dict: updated binding dict
        """
        if bindings[(edge[0], edge[1])] == 2:
            bindings[(edge[0], edge[1])] = 1
        else:
            bindings[(edge[0], edge[1])] = 2
        return bindings


    def restore(self, path: list, bindings: dict, bindings_list: list):
        """Helping function to restore all variables that have been changed when a node was added to the current path

        Args:
            path (list): current found path
            bindings (dict): current binding types of all edges in the graph
            bindings_list (list): current list of binding types (for each timestep) 

        Returns:
            int, list, dict, int, list: akt, path, bindings, bind, binding_list
        """
        akt = path[len(path)-2] # the node before the node that was wrongly added
        wrong = path[len(path)-1] # the node that was wrongly added
        path.pop(len(path)-1) # removing wrongly added node
        bindings = self.change_bindings(bindings, tuple(sorted((akt, wrong)))) # restore binding type of the edge that was wrongly traversed
        if len(path) == 1: # if only one node is in current path than there is no binding type of lastly traversed edge
            bind = 0
        else:
            if bindings[tuple(sorted((path[len(path)-2], path[len(path)-1])))] == 2: # restore bind
                bind = 1
            else:
                bind = 2
        bindings_list.pop(len(bindings_list)-1) # remove wrongly added binding dict
        return akt, path, bindings, bind, bindings_list


    def find_first_path(self, graph: nx.Graph, bindings: dict, end: int, path: list, akt: int, bind: int, bindings_list: list):
        """Find first soliton path between two exterior nodes
            A path can only be a soliton path if the edges traversed by the soliton have alternating binding types (1,2,1,2,..)

        Args:
            graph (nx.Graph): graph the path should be found in
            bindings (dict): current binding types of all edges in the graph
            end (int): end node of path
            path (list): current found path
            akt (int): node that was currently added to path
            bind (int): binding type of the last edge that was traversed
            bindings_list (list): current list of binding types (for each timestep)

        Returns:
            bool: if a path was found or not
            if a path was found, also:
                list: found path
                list: bindings_list for the found path
        """
        # base case: if end node is reachable then add end node to path and we are done
        if end in list(nx.neighbors(graph, akt)) and bindings[tuple(sorted((akt, end)))] != bind:
            path.append(end)
            bindings = self.change_bindings(bindings, (akt, end)) # change binding of traversed edge
            bindings_copy = bindings.copy() # make a copy of binding dict because otherwise we are using wrong references
            bindings_list.append(bindings_copy) # add binding dict of this timestep to bindings_list
            return True, path, bindings_list

        # iterate over all nodes that are adjacent to latest node in path
        for node in list(nx.neighbors(graph, akt)):
            if node != path[len(path)-2] and bindings[tuple(sorted((akt, node)))] != bind: # soliton is not allowed to make a direct turnaround and edge to next node has to have the right binding type
                path.append(node)
                bind = bindings[tuple(sorted((akt, node)))] # change bind to binding type of edge that was just traversed
                bindings = self.change_bindings(bindings, tuple(sorted((akt, node))))
                akt = node
                bindings_copy = bindings.copy()
                bindings_list.append(bindings_copy)
                # call function recursively: if we can find a path if we go further with the decision we just made (with the node we just added) then we are done
                if self.find_first_path(graph, bindings, end, path, akt, bind, bindings_list):
                    return True, path, bindings_list 
                else:
                    akt, path, bindings, bind, bindings_list = self.restore(path, bindings, bindings_list) # otherwise our decision was wrong, so we have to make it undone and try again with another node

        return False # if at some point no new node could be added, then no path can be found 
    

    def call_find_first_path(self):
        """Initialising some parameteres and then calling the find_first_path function with them

        Returns:
            bool: return False if no path could be found
            if a path could be found:
                list: found path
                list: bindings_list for the found path 
        """
        graph = self.soliton_graph.graph
        bindings = self.soliton_graph.bindings
        path = [self.start]
        akt = self.start
        bind = 0
        bindings_list = []
        bindings_copy = bindings.copy()
        bindings_list.append(bindings_copy)
        res = self.find_first_path(graph, bindings, self.end, path, akt, bind, bindings_list)
        if res != False:
            return res[1], res[2] # if a path was found, return the path and the bindings_list
        return res 


if __name__ == "__main__":

    my_graph = SolitonGraph('C1{1}=C{3}C1{=2}')
    #my_graph = SolitonGraph('C1=CC=CC=C1C{1}=CC{2}=CC=CC2=CC=CC=C2')
    automata = SolitonAutomata(my_graph, 1, 5)
    print(automata.soliton_path.adjacency_matrices_list)
    print(automata.soliton_path.path)
    #print(automata.call_find_first_path(7, 10))
