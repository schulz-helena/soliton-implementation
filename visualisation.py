"""Visualisation of a molecule as a graph
"""
import re

import matplotlib.pyplot as plt
import networkx as nx
from pysmiles import read_smiles
from rdkit import Chem
from rdkit.Chem import AllChem


def get_ring_systems(mol: Chem.MolFromSmiles):
    """Get all rings in a molecule

    Args:
        mol (Chem.MolFromSmiles): molecule

    Returns:
        list: all found rings (each ring in a seperate list)
    """
    rings = mol.GetRingInfo()
    systems = []
    for ring in rings.AtomRings():
        ring_ats = []
        for atom in ring:
            ring_ats.append(atom)
        systems.append(ring_ats)
    return systems


def get_ring_for_atom(rings: list, atom_id: int):
    """Get the ring the atom that was handed over is contained in

    Args:
        rings (list): list of all rings in molecule
        atom_id (int): id of the handed over atom

    Returns:
        list: ring that was found (None if no ring was found)
    """
    for ring in rings:
        if atom_id in ring:
            return ring


def double_edge(graph, node1, node2):
    pos = nx.get_node_attributes(graph, 'pos')
    coord1 = pos[node1]
    coord2 = pos[node2]
    x_values = [coord1[0], coord2[0]]
    y_values = [coord1[1], coord2[1]]

    dx = x_values[1]- x_values[0]
    dy = y_values[1]- y_values[0]
    slope = dy/dx

    #0.1 (the value that is added) is just some value for now, could be chosen better maybe 
    # TODO: compute slope for different edges and then find the perfect distance for them → find an algorithm/ formula to compute perfect distance
    # TODO: sometimes use + and sometimes - the distance (so second edge is always on outside or always on inside of circel)
    # TODO: make line of second edge thinner (so it has the same width as the normal edges in the graph)
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
        input_with_nums = re.sub(r"[C]", str(current), input_with_nums, count = 1)
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
                plt.plot(x_values, y_values, color = 'blue')
        else:
            if (bindings[(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())] == 2):
                graph.edges[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]['weight'] = 2
                graph.nodes[bond.GetBeginAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetBeginAtomIdx()] + 1
                graph.nodes[bond.GetEndAtomIdx()]['weight'] = nx.get_node_attributes(graph, 'weight')[bond.GetEndAtomIdx()] + 1
                x_values, y_values = double_edge(graph, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                plt.plot(x_values, y_values, color = 'blue')
    return graph


def mol_to_nx(mol: Chem.MolFromSmiles):
    """Convert rdkit object (molecule) to networkx object (graph)

    Args:
        mol (Chem.MolFromSmiles): molecule

    Returns:
        nx.DiGraph: graph
    """
    graph = nx.DiGraph()
    AllChem.Compute2DCoords(mol)
    rings = get_ring_systems(mol)
    arom_atoms = set()
    node_label = 'a'

    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        x_coord = pos.x
        y_coord = pos.y * 1.2
        graph.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   num=node_label,
                   pos=(x_coord, y_coord))
        node_label = chr(ord(node_label)+1)  # noch ändern: nach z dann aa usw.

    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
        # if its a double bond add a second edge (in other direction)
        if (bond.GetBondType() == Chem.BondType.DOUBLE):
            graph.add_edge(bond.GetEndAtomIdx(),
                       bond.GetBeginAtomIdx(),
                       bond_type=bond.GetBondType())
        # in aromatic rings every second edge is a double bond
        # if it is aromatic get corresponding ring and add every second atom from that ring to set (arom_atoms)
        if (bond.GetBondType() == Chem.BondType.AROMATIC):
            ring = get_ring_for_atom(rings, bond.GetBeginAtomIdx())
            if (ring is not None):
                for i in range(0, len(ring), 2):
                    arom_atoms.add(ring[i])
                rings.remove(ring)
        # if beginning atom of bond is in arom_atoms add a second edge (in other direction)
        if (bond.GetBeginAtomIdx() in arom_atoms):
            graph.add_edge(bond.GetEndAtomIdx(),
                       bond.GetBeginAtomIdx(),
                       bond_type=bond.GetBondType())
    return graph

def mol_to_nx2(mol: Chem.MolFromSmiles):
    """Convert rdkit object (molecule) to networkx object (graph)

    Args:
        mol (Chem.MolFromSmiles): molecule

    Returns:
        nx.Graph: graph
    """
    graph = nx.Graph()
    AllChem.Compute2DCoords(mol)
    rings = get_ring_systems(mol)
    arom_atoms = set()
    node_label = 'a'

    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        x_coord = pos.x
        y_coord = pos.y * 1.2
        graph.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   num=node_label,
                   pos=(x_coord, y_coord))
        node_label = chr(ord(node_label)+1)  # noch ändern: nach z dann aa usw.

    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
        # if its a double bond add a second edge (in other direction)
        if (bond.GetBondType() == Chem.BondType.DOUBLE):
            x_values, y_values = double_edge(graph, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            plt.plot(x_values, y_values, color = 'blue')
        # in aromatic rings every second edge is a double bond
        # if it is aromatic get corresponding ring and add every second atom from that ring to set (arom_atoms)
        if (bond.GetBondType() == Chem.BondType.AROMATIC):
            ring = get_ring_for_atom(rings, bond.GetBeginAtomIdx())
            if (ring is not None):
                for i in range(0, len(ring), 2):
                    arom_atoms.add(ring[i])
                rings.remove(ring)
        # if beginning atom of bond is in arom_atoms add a second edge (in other direction)
        if (bond.GetBeginAtomIdx() in arom_atoms):
            x_values, y_values = double_edge(graph, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            plt.plot(x_values, y_values, color = 'blue')
    return graph

if __name__ == "__main__":

    #SMILES = 'C1=CC=CC=C1C=CC=CC2=CC=CC=C2'
    #SMILES = 'C1=C(C)C=CC=C1C(C)=C(C)C=CC2=CC(C)=C(C)C=C2' #mit externen Knoten
    #SMILES = 'C1=CC=C2C=CC=CC2=C1' #zwei Ringe aneinander
    #mol = Chem.MolFromSmiles(SMILES)

    #molecule_nx = mol_to_nx3('C1=CC=CC=C1C=CC=CC2=CC=CC=C2')
    molecule_nx = mol_to_nx3('C1=CC=CC=C1C=CC=C{=1}C2=CC=CC=C2')

    #molecule_atom = nx.get_node_attributes(molecule_nx, 'num')
    labels = nx.get_node_attributes(molecule_nx, 'label')
    pos = nx.get_node_attributes(molecule_nx, 'pos')

    plt.axis('equal')

    #x_values, y_values = double_edge(molecule_nx, 1, 2)
    #plt.plot(x_values, y_values, color = 'blue')

    '''nx.draw(molecule_nx,
            pos,
            labels=molecule_atom,
            with_labels=True,
            node_color='white',
            connectionstyle='bar, fraction=0.03',
            arrowstyle='-')'''

    print(nx.get_node_attributes(molecule_nx, 'weight'))
    print(nx.get_edge_attributes(molecule_nx, 'weight'))

    nx.draw(molecule_nx,
            pos,
            labels=labels,
            with_labels=True,
            node_color='white')

    plt.show()

    # adjacency matrix
    '''matrix = nx.to_numpy_matrix(molecule_nx)
    print(matrix)'''
