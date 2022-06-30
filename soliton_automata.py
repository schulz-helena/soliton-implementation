"""Class that represents a soliton automata and therefore finds all soliton paths for a given pair of exterior nodes 

    Attributes: soliton graph the automata is based on, start and end node for paths, list of soliton paths (right now only the first found soliton path)
"""
import copy

import networkx as nx

from soliton_graph import SolitonGraph
from soliton_path import SolitonPath


class SolitonAutomata:

    def __init__(self, soliton_graph, start, end):
        #path, bindings_list = self.call_find_first_path()
        #self.soliton_path = SolitonPath(path, self.soliton_graph, bindings_list)
        self.soliton_graph = soliton_graph
        self.start = self.soliton_graph.exterior_nodes_reverse[str(start)]
        self.end = self.soliton_graph.exterior_nodes_reverse[str(end)]
        self.paths_ids = self.call_find_all_paths()
        #if self.paths_ids == []:
            #print("There exists no soliton path between these exterior nodes")
        self.paths = [] # paths with node labels instead of node ids
        self.paths_for_user = [] # representation of the path the user gets as an output
        for path_ids in self.paths_ids:
            path = copy.deepcopy(path_ids)
            path_string = "" # for user output
            for i in range(0, len(path)):
                path[i] = nx.get_node_attributes(self.soliton_graph.graph, 'label')[path[i]]
                if i == len(path) - 1:
                    path_string = path_string + f"{path[i]}"
                else:
                    path_string = path_string + f"{path[i]} - "
            self.paths.append(path)
            self.paths_for_user.append(path_string)


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


    def build_copies(self, akt: int, path: list, bindings: dict, bind: int):
        """Helping function to copy all variables that are changed during find_all_paths

        Args:
            akt (int): node that was currently added to path
            path (list): current found path
            bindings (dict): current binding types of all edges in the graph
            bind (int): binding type of the last edge that was traversed

        Returns:
            int, list, dict, int: akt_copy, path_copy, bindings_copy, bind_copy
        """
        akt_copy = copy.deepcopy(akt)
        path_copy = copy.deepcopy(path)
        bindings_copy = copy.deepcopy(bindings)
        bind_copy = copy.deepcopy(bind)
        return akt_copy, path_copy, bindings_copy, bind_copy


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


    def find_all_paths(self, graph: nx.Graph, bindings: dict, end: int, path: list, akt: int, bind: int, paths: list):
        """Find all possible soliton paths between two exterior nodes
            A path can only be a soliton path if the edges traversed by the soliton have alternating binding types (1,2,1,2,..)

        Args:
            graph (nx.Graph): graph the path should be found in
            bindings (dict): current binding types of all edges in the graph
            end (int): end node of path
            path (list): current found path
            akt (int): node that was currently added to path
            bind (int): binding type of the last edge that was traversed
            paths (list): contains all currently found paths

        Returns:
            list: contains all found paths (is empty if no path exists)
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
    

    def call_find_first_path(self):
        """Initialising some parameteres and then calling the find_first_path function with them

        Returns:
            bool: return False if no path could be found
            if a path could be found:
                list: found path
                list: bindings_list for the found path 
        """
        soliton_graph_copy = copy.deepcopy(self.soliton_graph) # working on copy of graph so no unwanted changes are made
        graph = soliton_graph_copy.graph
        bindings = soliton_graph_copy.bindings
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


    def call_find_all_paths(self):
        """Initialising some parameteres and then calling the find_all_paths function with them

        Returns:
            list: list of all found paths (returns empty list when no path is found)
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



if __name__ == "__main__":

    #my_graph = SolitonGraph('C1{1}=C{3}C1{=2}')
    #my_graph = SolitonGraph('C1=CC=CC=C1C=C{1}C=CC{2}=CC2=CC=CC=C2')
    my_graph = SolitonGraph('C1=CC{1}=CC=CC=C1C{2}=CC=CC=CC=CC=CC=C{3}C2=C{4}C=CC=CC=C2')
    print(my_graph.labels)
    print(my_graph.exterior_nodes)
    automata = SolitonAutomata(my_graph, 2, 3)
    #automata = SolitonAutomata(my_graph, 1, 1)
    #print(automata.soliton_path.adjacency_matrices_list)
    #print(automata.paths_ids)
    print(automata.paths)
