import re

import matplotlib.pyplot as plt
import networkx as nx
from pysmiles import read_smiles
from rdkit import Chem
from rdkit.Chem import AllChem


def user_input_to_dict(user_input: str):
    # unfinished and unused
    iterable_string = user_input
    dict_nodes = {}
    dict_edges = {}
    current = 0
    node_label = 'a'
    while True:
        if iterable_string == "":
            break
        element = re.search(r"[C][)]*[0-9]*[(]*[=-]*", iterable_string).group()
        dict_nodes[current] = node_label
        # increase counter variables and delete found chars from string
        node_label = chr(ord(node_label)+1)
        current += 1
        iterable_string = re.sub(r"[C][)]*[0-9]*[(]*[=-]*", "", iterable_string, count = 1)
        print(element)
        print(iterable_string)
        print(dict_nodes)

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

#using a directed graph
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

#using an undirected graph but imprecise way of finding/ drawing double edges
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

#Update function with blue colored nodes instead of black pebble as soliton on a node
def simple_update(frame_num: int, graph: nx.Graph, ax: Axes, path: list, bindings: dict, xs_and_ys: dict):
    """Update function that changes graph while soliton animation

    Args:
        frame_num (int): frame number (is increased every time animation calls this function)
        G (nx.Graph): graph that is supposed to be animated
        ax (Axes): axes of the plot
        path (list): path that the soliton should traverse
    """
    ax.clear()
    ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

    label = nx.get_node_attributes(graph, 'label')
    pos = nx.get_node_attributes(graph, 'pos')
    color_map = []
    for node in graph:
        if node != path[frame_num]:
            color_map.append('white')
        else:
            color_map.append('blue')
            if frame_num != 0:
                if (path[frame_num-1], path[frame_num]) in bindings:
                    if bindings[(path[frame_num-1], path[frame_num])] == 2:
                        bindings[(path[frame_num-1], path[frame_num])] = 1
                    else:
                        bindings[(path[frame_num-1], path[frame_num])] = 2
                else:
                    if bindings[(path[frame_num], path[frame_num-1])] == 2:
                        bindings[(path[frame_num], path[frame_num-1])] = 1
                    else:
                        bindings[(path[frame_num], path[frame_num-1])] = 2
    
    for edge in graph.edges:
        if bindings[edge] == 2:
            x_values, y_values = xs_and_ys[edge]
            plt.plot(x_values, y_values, color = 'black', linewidth = 1)

    nx.draw(graph, pos=pos, labels=label, with_labels=True, node_color = color_map, ax=ax)

    # Set the title
    ax.set_title(f"Frame {frame_num}")

if __name__ == "__main__":

    SMILES = 'C1=CC=CC=C1C=CC=CC2=CC=CC=C2'
    #SMILES = 'C1=CC=C2C=CC=CC2=C1' #zwei Ringe aneinander
    mol = Chem.MolFromSmiles(SMILES)

    molecule_nx = mol_to_nx2('C1=CC=CC=C1C=CC=CC2=CC=CC=C2')

    molecule_atom = nx.get_node_attributes(molecule_nx, 'num')
    pos = nx.get_node_attributes(molecule_nx, 'pos')

    plt.axis('equal')

    '''nx.draw(molecule_nx,
            pos,
            labels=molecule_atom,
            with_labels=True,
            node_color='white',
            connectionstyle='bar, fraction=0.03',
            arrowstyle='-')'''

    nx.draw(molecule_nx,
            pos,
            labels=molecule_atom,
            with_labels=True,
            node_color='white')

    plt.show()

    # adjacency matrix
    '''matrix = nx.to_numpy_matrix(molecule_nx)
    print(matrix)'''
