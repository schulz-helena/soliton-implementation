"""Class that represents a soliton graph

    Attributes: exterior nodes (node id as key, node label as value and other way around), two smiles strings (for use with pysmiles and for use with rdkit), 
                binding types (edge as key, binding type as value), nx-graph, positions for a second line for each edge 
"""
import re

import networkx as nx
from pysmiles import read_smiles
from rdkit import Chem
from rdkit.Chem import AllChem


class SolitonGraph:

    def __init__(self, user_input: str):
        self.exterior_nodes, self.exterior_nodes_reverse = self.find_exterior_nodes(user_input)
        self.pysmiles_smiles = self.create_pysmiles_smiles(user_input)
        self.rdkit_smiles = self.create_rdkit_smiles(user_input)
        self.bindings = self.create_binding_dict()
        self.graph = self.smiles_to_graph()
        self.double_edge_positions = self.find_double_edge_positions()
        self.labels = nx.get_node_attributes(self.graph, 'label')


    def set_bindings(self, bindings):
        self.bindings = bindings
        self.graph = self.smiles_to_graph()


    def find_exterior_nodes(self, user_input: str):
        """Transform the users input into dictionary with exterior nodes

        Args:
            user_input (str): user input

        Returns:
            dict: exterior nodes with node ids as keys as node labels as values
            dict: exterior nodes with node labels as keys as node ids as values
        """
        # exterior nodes are put in "{}" in user input
        exterior_nodes = {} # dictionary for exterior nodes
        exterior_nodes_reverse = {}
        current = 0
        # find node labels of exterior nodes (numbers)
        matches_labels = re.findall(r"[{][-=]*[0-9]*[}]", user_input)
        # replace node labels with Cs, then count Cs in string and replace each C with count
        input_with_c = re.sub(r"[{][-=]*[0-9]*[}]", "{C}", user_input)
        input_with_nums = input_with_c
        while True:
            if (re.search(r"[C]", input_with_nums) is None):
                break
            input_with_nums = re.sub(r"[CSNOF]", str(current), input_with_nums, count = 1) # [CSNOF]
            current += 1
        # find node ids in input_with_nums (because node id is just the atom count)
        matches_ids = re.findall(r"[{][-=]*[0-9]*[}]", input_with_nums)
        # put exterior nodes in dictionary (node id as key and node label as value)
        for i in range (0, len(matches_labels)):
            matches_ids[i] = re.sub(r"[}]", "", re.sub(r"[{]", "", matches_ids[i]))
            matches_labels[i] = re.sub(r"[}]", "", re.sub(r"[{][-=]*", "", matches_labels[i]))
            exterior_nodes[int(matches_ids[i])] = matches_labels[i]
            exterior_nodes_reverse[matches_labels[i]] = int(matches_ids[i])

        return exterior_nodes, exterior_nodes_reverse


    def create_pysmiles_smiles(self, user_input: str):
        """Transform user input in smiles representation (treating exterior nodes as Cs now)

        Args:
            user_input (str): user input

        Returns:
            str: smiles string (used with pysmiles)
        """
        pysmiles_smiles = re.sub(r"[{][0-9]*[}]", "(C)", user_input)
        pysmiles_smiles = re.sub(r"[{][-][0-9]*[}]", "(-C)", pysmiles_smiles)
        pysmiles_smiles = re.sub(r"[{][=][0-9]*[}]", "(=C)", pysmiles_smiles)
        return pysmiles_smiles


    def create_rdkit_smiles(self, user_input: str):
        """Transform user input in extra smiles representation
            because rdkit needs smiles without double edges at exterior nodes (otherwise some valence error occurs)

        Args:
            user_input (str): user input

        Returns:
            str: modified smiles string (used with rdkit)
        """
        rdkit_smiles = re.sub(r"[{][-=]*[0-9]*[}]", "(C)", user_input)
        return rdkit_smiles


    def create_binding_dict(self):
        """Build a dictionary that contains the binding type for each edge (1 for single, 2 for double binding)

        Returns:
            dict: bindings
        """
        mol_pysmiles = read_smiles(self.pysmiles_smiles, reinterpret_aromatic=False) # binding information are taken from pysmiles
        bindings = nx.get_edge_attributes(mol_pysmiles, 'order')
        bindings_sorted_tuples = {}
        for edge in bindings:
            val = bindings[edge]
            bindings_sorted_tuples[tuple(sorted(edge))] = val
        return bindings_sorted_tuples


    def next_node_label(self, node_label: str):
        """Helping function to find the next node label for a given node label (for initialisation of grapg)

        Args:
            node_label (str): given node label

        Returns:
            str: the computed (next) node label
        """
        if len(node_label) == 1: # e.g.: a -> b
            node_label = chr(ord(node_label)+1)
        elif node_label[1] == 'z': # e.g. bz -> ca
            node_label_list = list(node_label) # convert to list so we can change chars at certain index
            node_label_list[0] = chr(ord(node_label[0])+1)
            node_label_list[1] = 'a'
            node_label = "".join(node_label_list) # convert back to string
        else: #e.g. ab -> ac
            node_label_list = list(node_label)
            node_label_list[1] = chr(ord(node_label[1])+1)
            node_label = "".join(node_label_list)
        return node_label


    def smiles_to_graph(self):
        """Transform user input into molecule and then into nx graph

        Returns:
            nx.Graph: graph that visualizes the molecule
        """
        mol_rdkit = Chem.MolFromSmiles(self.rdkit_smiles) # atom position information are taken from rdkit
        AllChem.Compute2DCoords(mol_rdkit)
        graph = nx.Graph()
        if (len(mol_rdkit.GetAtoms()) - len(self.exterior_nodes)) > 26: # if we have more than 26 atoms then node labels a - z are not sufficient
            node_label = 'aa'
        else:
            node_label = 'a'

        for atom in mol_rdkit.GetAtoms():
            pos = mol_rdkit.GetConformer().GetAtomPosition(atom.GetIdx())
            x_coord = pos.x
            y_coord = pos.y * 1.2
            if atom.GetIdx() in self.exterior_nodes:
                graph.add_node(atom.GetIdx(),
                    label=self.exterior_nodes[atom.GetIdx()],
                    pos=(x_coord, y_coord),
                    weight = 0)
            else:
                graph.add_node(atom.GetIdx(),
                    label=node_label,
                    pos=(x_coord, y_coord),
                    weight = 0)
                node_label = self.next_node_label(node_label)

        for bond in mol_rdkit.GetBonds():
            graph.add_edge(bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                    weight = 1)
            # increase node weight of begin atom and end atom by 1
            graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
            graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
            # if its a double bond increase edge and node weights
            if ((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in self.bindings):
                if (self.bindings[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] == 2):
                    # put edge weight to 2
                    graph.edges[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]['weight'] = 2
                    # increase node weight of begin atom and end atom by 1
                    graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
                    graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
            elif ((bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) in self.bindings):
                if (self.bindings[(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())] == 2):
                    graph.edges[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]['weight'] = 2
                    graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
                    graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
        return graph


    def find_double_edge_positions(self):
        """Build a dictionary that contains the positions for a second line for each edge

        Returns:
            dict: double edge positions
        """

        def double_edge(graph: nx.Graph, node1: int, node2: int):
            """Get x and y values for drawing a line as a second edge

            Args:
                graph (nx.Graph): graph the second edge should be added to
                node1 (int): beginning node of edge
                node2 (int): end node of edge

            Returns:
                list: x values
                list: y values
            """
            pos = nx.get_node_attributes(graph, 'pos')
            coord1 = pos[node1]
            coord2 = pos[node2]
            x_values = [coord1[0], coord2[0]]
            y_values = [coord1[1], coord2[1]]

            dx = x_values[1] - x_values[0]
            # no division by 0 allowed
            if (dx == 0):
                dx = 0.001
            dy = y_values[1] - y_values[0]
            slope = dy/dx

            # 0.1 (the value that is added) is just some value for now, could be chosen better maybe 
            # TODO: compute slope for different edges and then find the perfect distance for them → find an algorithm/ formula to compute perfect distance
            # TODO: sometimes use + and sometimes - the distance (so second edge is always on outside or always on inside of circel)
            if (abs(slope) >= 1):
                x_values = [coord1[0]+0.1, coord2[0]+0.1]
                y_values = [coord1[1], coord2[1]]
            else:
                x_values = [coord1[0], coord2[0]]
                y_values = [coord1[1]+0.1, coord2[1]+0.1]

            return x_values, y_values

        double_edge_positions = {}
        for edge in self.graph.edges:
            x_values, y_values = double_edge(self.graph, edge[0], edge[1])
            double_edge_positions[edge] = x_values, y_values
        return double_edge_positions
    

    def validate_soliton_graph(self):
        """Check if graph is a soliton graph
        """
        errors = []
        weights = nx.get_node_attributes(self.graph, 'weight')
        #labels = nx.get_node_attributes(self.graph, 'label')
        # No self-loops
        selfloops = list(nx.nodes_with_selfloops(self.graph))
        if len(selfloops) > 0:
            for node in selfloops:
                errors.append(f"Self-loop at node {self.labels[node]}")
                #print(f"Self-loop at node {labels[node]}")
        # Only node degress between 1 and 3 allowed
        for (node, val) in self.graph.degree():
            if val > 3:
                errors.append(f"Node {self.labels[node]} has too many neighbours")
                #print(f"Node {labels[node]} has too many neighbours")
        # exterior nodes must have weight of 1 or 2 and must have degree 1
        for key in self.exterior_nodes:
            if weights[key] > 2:
                errors.append(f"The weight of node {self.labels[key]} is too high")
                #print(f"The weight of node {labels[key]} is too high")
            del weights[key]
            if self.graph.degree(key) > 1:
                errors.append(f"Node {self.labels[key]} has too many neighbours")
                #print(f"Node {labels[key]} has too many neighbours")
        # Inner nodes have exactly one double edge
        for node in weights:
            if weights[node] > self.graph.degree(node) + 1:
                errors.append(f"The weight of node {self.labels[node]} is too high")
                #print(f"The weight of node {labels[node]} is too high")
            elif weights[node] < self.graph.degree(node) + 1:
                errors.append(f"The weight of node {self.labels[node]} is too low")
                #print(f"The weight of node {labels[node]} is too low")
        # There has to be at least one exterior node
        if len(self.exterior_nodes) < 1:
            errors.append("You must have at least one exterior node")
            #print("You must have at least one exterior node")
        return errors
