"""Visualisation of a molecule as a graph
"""
import re

import matplotlib.pyplot as plt
import networkx as nx
from pysmiles import read_smiles
from rdkit import Chem
from rdkit.Chem import AllChem


def double_edge(graph, node1, node2):
    pos = nx.get_node_attributes(graph, 'pos')
    coord1 = pos[node1]
    coord2 = pos[node2]
    x_values = [coord1[0], coord2[0]]
    y_values = [coord1[1], coord2[1]]

    dx = x_values[1] - x_values[0]
    #no division by 0 allowed
    if (dx == 0):
        dx = 0.001
    dy = y_values[1] - y_values[0]
    slope = dy/dx

    #0.1 (the value that is added) is just some value for now, could be chosen better maybe 
    # TODO: compute slope for different edges and then find the perfect distance for them → find an algorithm/ formula to compute perfect distance
    # TODO: sometimes use + and sometimes - the distance (so second edge is always on outside or always on inside of circel)
    if (abs(slope) >= 1):
        x_values = [coord1[0]+0.1, coord2[0]+0.1]
        y_values = [coord1[1], coord2[1]]
    else:
        x_values = [coord1[0], coord2[0]]
        y_values = [coord1[1]+0.1, coord2[1]+0.1]

    return x_values, y_values


def transform_user_input(user_input: str):
    #external nodes are put in "{}" in user input
    external_nodes = {} #dictionary for external nodes
    current = 0
    #find node labels of external nodes (numbers)
    matches_labels = re.findall(r"[{][-=]*[0-9]*[}]", user_input)
    #replace node labels with Cs, then count Cs in string and replace each C with count
    input_with_C = re.sub(r"[{][-=]*[0-9]*[}]", "{C}", user_input)
    input_with_nums = input_with_C
    while True:
        if (re.search(r"[C]", input_with_nums) == None):
            break
        input_with_nums = re.sub(r"[CSNOF]", str(current), input_with_nums, count = 1) #TODO: also put in other possible atoms
        current += 1
    #find node ids in input_with_nums (because node id is just the atom count)
    matches_ids = re.findall(r"[{][-=]*[0-9]*[}]", input_with_nums)
    #put external nodes in dictionary (node id as key and node label as value)
    for i in range (0, len(matches_labels)):
        matches_ids[i] = re.sub(r"[}]", "", re.sub(r"[{]", "", matches_ids[i]))
        matches_labels[i] = re.sub(r"[}]", "", re.sub(r"[{][-=]*", "", matches_labels[i]))
        external_nodes[int(matches_ids[i])] = int(matches_labels[i])
    #transform user input in smiles representation (treating external nodes as Cs now)
    input_as_smiles = re.sub(r"[{][0-9]*[}]", "(C)", user_input)
    input_as_smiles = re.sub(r"[{][-][0-9]*[}]", "(-C)", input_as_smiles)
    input_as_smiles = re.sub(r"[{][=][0-9]*[}]", "(=C)", input_as_smiles)

    #rdkit needs smiles without double edges at external nodes (otherwise some valence error occurs)
    #but binding information are taken from pysmiles anyway
    rdkit_smiles = re.sub(r"[{][-=]*[0-9]*[}]", "(C)", user_input)

    return external_nodes, input_as_smiles, rdkit_smiles


def mol_to_nx3(user_input: str):
    external_nodes, smiles, rdkit_smiles = transform_user_input(user_input)
    mol_pysmiles = read_smiles(smiles, reinterpret_aromatic=False)
    bindings = nx.get_edge_attributes(mol_pysmiles, 'order')
    mol_rdkit = Chem.MolFromSmiles(rdkit_smiles)
    AllChem.Compute2DCoords(mol_rdkit)
    graph = nx.Graph()
    node_label = 'a'

    for atom in mol_rdkit.GetAtoms():
        pos = mol_rdkit.GetConformer().GetAtomPosition(atom.GetIdx())
        x_coord = pos.x
        y_coord = pos.y * 1.2
        if atom.GetIdx() in external_nodes:
            graph.add_node(atom.GetIdx(),
                   label=str(external_nodes[atom.GetIdx()]),
                   pos=(x_coord, y_coord),
                   weight = 0)
        else:
            graph.add_node(atom.GetIdx(),
                   label=node_label,
                   pos=(x_coord, y_coord),
                   weight = 0)
            node_label = chr(ord(node_label)+1)  # TODO: noch ändern: nach z dann aa usw.

    for bond in mol_rdkit.GetBonds():
        #print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        graph.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   weight = 1)
        #increase node weight of begin atom and end atom by 1
        graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
        graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
        # if its a double bond add a second edge (in other direction)
        if ((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in bindings):
            if (bindings[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] == 2):
                #put edge weight to 2
                graph.edges[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]['weight'] = 2
                graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
                graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
                x_values, y_values = double_edge(graph, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                plt.plot(x_values, y_values, color = 'black', linewidth = 1)
        elif ((bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) in bindings):
            if (bindings[(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())] == 2):
                graph.edges[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]['weight'] = 2
                graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
                graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
                x_values, y_values = double_edge(graph, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                plt.plot(x_values, y_values, color = 'black', linewidth = 1)
    return graph


if __name__ == "__main__":

    molecule_nx = mol_to_nx3('COC1=CC2=C(C=C1)N(C(=O)C1=CC=C(Cl)C=C1)C(C)=C2CC(=O)OCC(O)=O')

    labels = nx.get_node_attributes(molecule_nx, 'label')
    pos = nx.get_node_attributes(molecule_nx, 'pos')

    plt.axis('equal')

    #print(nx.get_node_attributes(molecule_nx, 'weight'))
    #print(nx.get_edge_attributes(molecule_nx, 'weight'))
    # adjacency matrix
    #matrix = nx.to_numpy_matrix(molecule_nx)
    #print(matrix)

    nx.draw(molecule_nx,
            pos,
            labels=labels,
            with_labels=True,
            node_color='white')

    plt.show()
