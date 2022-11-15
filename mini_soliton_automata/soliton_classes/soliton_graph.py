"""Soliton graph.
"""
import re

import networkx as nx
from pysmiles import read_smiles
from rdkit import Chem
from rdkit.Chem import AllChem


class SolitonGraph:
    """Representation of a soliton graph.
    """

    def __init__(self, user_input: str):
        """Initializes a `SolitonGraph` object by using the input string of the user.
        """
        self.exterior_nodes: dict
        """Exterior nodes (node id as key, node label as value)."""
        self.exterior_nodes_reverse: dict
        """Exterior nodes (node label as key, node id as value)."""
        self.exterior_nodes, self.exterior_nodes_reverse = self.find_exterior_nodes(user_input)
        self.pysmiles_smiles: str = self.create_pysmiles_smiles(user_input)
        """`SMILES` string for use with `pysmiles`."""
        self.rdkit_smiles: str = self.create_rdkit_smiles(user_input)
        """`SMILES` string for use with `rdkit`."""
        self.bindings: dict = self.create_binding_dict()
        """Binding types (edge as key, binding type as value)."""
        self.graph: nx.Graph = self.smiles_to_graph()
        """Graph that represents the molecule."""
        self.double_edge_positions: dict = self.find_double_edge_positions()
        """Positions for a second line (that can be plotted) for each edge."""
        self.labels: dict = nx.get_node_attributes(self.graph, 'label')
        """Node labels (node id as key, node label as value)."""
        self.way: str = user_input
        """Way that led to the soliton graph (consists of user input of initial soliton graph of a soliton automata and soliton paths)."""


    def set_bindings(self, bindings: dict):
        """Sets the soliton graph's bindings to a new binding dictionary.

        Args:
            bindings (dict): New binding dictionary.
        """
        self.bindings = bindings
        self.graph = self.smiles_to_graph()


    def find_exterior_nodes(self, user_input: str):
        """Transforms the users input into dictionary with exterior nodes.

        Args:
            user_input (str): User input.

        Returns:
            dict: Exterior nodes with node ids as keys and node labels as values.
            dict: Exterior nodes with node labels as keys and node ids as values.
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
            input_with_nums = re.sub(r"[C]", str(current), input_with_nums, count = 1)
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

    def exterior_nodes_name_collision(self):
        """Checks for name collisions between exterior nodes.

        Returns:
            bool: `True`, if user used same exterior node label more than once, `False` otherwise.
        """
        flipped = {}
        for key, value in self.exterior_nodes.items():
            if value not in flipped:
                flipped[value] = [key]
            else:
                return True

        return False

    def create_pysmiles_smiles(self, user_input: str):
        """Transforms user input in `SMILES` representation (treating exterior nodes as Cs now).

        Args:
            user_input (str): User input.

        Returns:
            str: `SMILES` string (used with `pysmiles`).
        """
        pysmiles_smiles = re.sub(r"[{][0-9]*[}]", "(C)", user_input)
        pysmiles_smiles = re.sub(r"[{][-][0-9]*[}]", "(-C)", pysmiles_smiles)
        pysmiles_smiles = re.sub(r"[{][=][0-9]*[}]", "(=C)", pysmiles_smiles)

        return pysmiles_smiles


    def create_rdkit_smiles(self, user_input: str):
        """Transforms user input in extra `SMILES` representation
            because `rdkit` needs string without double edges at exterior nodes (otherwise some valence error occurs).

        Args:
            user_input (str): User input.

        Returns:
            str: Modified `SMILES` string (used with `rdkit`).
        """
        rdkit_smiles = re.sub(r"[{][-=]*[0-9]*[}]", "(C)", user_input)

        return rdkit_smiles


    def create_binding_dict(self):
        """Builds a dictionary that contains the binding type for each edge (1 for single, 2 for double binding).

        Returns:
            dict: Bindings (where the two nodes of the edge are sorted).
        """
        mol_pysmiles = read_smiles(self.pysmiles_smiles, reinterpret_aromatic=False) # binding information are taken from pysmiles (to ignore aromaticity)
        bindings = nx.get_edge_attributes(mol_pysmiles, 'order')
        bindings_sorted_tuples = {}
        for edge in bindings:
            val = bindings[edge]
            bindings_sorted_tuples[tuple(sorted(edge))] = val

        return bindings_sorted_tuples


    def next_node_label(self, node_label: str):
        """Finds the next node label for a given node label. Used in initialisation of graph in `smiles_to_graph`.

        Args:
            node_label (str): Given node label.

        Returns:
            str: Computed (next) node label.
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
        """Transforms user input into `rdkit` molecule and then into `networkx` graph.

        Returns:
            nx.Graph: Graph that visualizes the molecule.
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
        """Builds a dictionary that contains the positions for a second line for each edge.

        Returns:
            dict: Double edge positions.
        """

        def double_edge(graph: nx.Graph, node1: int, node2: int):
            """Gets x and y values for drawing a line as a second edge.

            Args:
                graph (nx.Graph): Graph the second edge should be added to.
                node1 (int): Beginning node of edge.
                node2 (int): End node of edge.

            Returns:
                list: X values.
                list: Y values.
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

            # distance between the two edges (0.11 as default, other values for extremely small or large graphs):
            distance = 0.11
            if len(pos) == 2:
                distance = 0.03
            elif len(pos) == 3:
                distance = 0.05
            elif len(pos) == 4:
                distance = 0.08
            elif len(pos) > 4 and len(pos) < 7:
                distance = 0.09
            elif len(pos) >= 26 and len(pos) < 42:
                distance = 0.15
            elif len(pos) >= 42:
                distance = 0.2

            if (abs(slope) >= 1):
                x_values = [coord1[0]+distance, coord2[0]+distance]
                y_values = [coord1[1], coord2[1]]
            else:
                x_values = [coord1[0], coord2[0]]
                y_values = [coord1[1]+distance, coord2[1]+distance]

            return x_values, y_values

        double_edge_positions = {}
        for edge in self.graph.edges:
            x_values, y_values = double_edge(self.graph, edge[0], edge[1])
            double_edge_positions[edge] = x_values, y_values

        return double_edge_positions
    

    def validate_soliton_graph(self):
        """Checks if graph is a soliton graph.

        Returns:
            list: All the problems that keep the graph from being a soliton graph.
        """
        errors = []
        weights = nx.get_node_attributes(self.graph, 'weight')
        # No self-loops
        selfloops = list(nx.nodes_with_selfloops(self.graph))
        if len(selfloops) > 0:
            for node in selfloops:
                errors.append(f"Self-loop at node {self.labels[node]}")
        # Only node degrees between 1 and 3 allowed
        for (node, val) in self.graph.degree():
            if val > 3:
                errors.append(f"Node {self.labels[node]} has too many neighbours")
        # exterior nodes must have weight of 1 or 2 and must have degree 1
        for key in self.exterior_nodes:
            if weights[key] > 2:
                errors.append(f"The weight of node {self.labels[key]} is too high")
            del weights[key]
            if self.graph.degree(key) > 1:
                errors.append(f"Node {self.labels[key]} has too many neighbours")
        # Inner nodes have exactly one double edge
        for node in weights:
            if weights[node] > self.graph.degree(node) + 1:
                errors.append(f"The weight of node {self.labels[node]} is too high")
            elif weights[node] < self.graph.degree(node) + 1:
                errors.append(f"The weight of node {self.labels[node]} is too low")
        # There has to be at least one exterior node
        if len(self.exterior_nodes) < 1:
            errors.append("You must have at least one exterior node")

        return errors
