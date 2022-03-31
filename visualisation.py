import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem


def GetRingSystems(mol):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = []
        for atom in ring:
            ringAts.append(atom)
        systems.append(ringAts)
    return systems


def getRingForAtom(rings, atomId):
    for ring in rings:
        if atomId in ring:
            return ring

# converting rdkit object to networkx object
def mol_to_nx(mol):
    G = nx.DiGraph()
    AllChem.Compute2DCoords(mol)
    rings = GetRingSystems(mol)
    aromAtoms = set()
    nodeLabel = 'a'

    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xCoord = pos.x
        yCoord = pos.y * 0.02
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   num=nodeLabel,
                   pos=(xCoord, yCoord))
        nodeLabel = chr(ord(nodeLabel)+1)  # noch ändern: nach z dann aa usw.

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
        # if its a double bond add a second edge (in other direction)
        if (bond.GetBondType() == Chem.BondType.DOUBLE):
            G.add_edge(bond.GetEndAtomIdx(),
                       bond.GetBeginAtomIdx(),
                       bond_type=bond.GetBondType())
        # in aromatic rings every second edge is a double bond
        # if it is aromatic get corresponding ring and add every second atom from that ring to set (aromAtoms)
        if (bond.GetBondType() == Chem.BondType.AROMATIC):
            ring = getRingForAtom(rings, bond.GetBeginAtomIdx())
            if (ring != None):
                for i in range(0, len(ring), 2):
                    aromAtoms.add(ring[i])
                rings.remove(ring)
        # if beginning atom of bond is in aromAtoms add a second edge (in other direction)
        if (bond.GetBeginAtomIdx() in aromAtoms):
            G.add_edge(bond.GetEndAtomIdx(),
                       bond.GetBeginAtomIdx(),
                       bond_type=bond.GetBondType())
    return G


if __name__ == "__main__":

    smiles = 'C1=CC=CC=C1C=CC=CC2=CC=CC=C2'
    #smiles = 'C1(=CC=CC=C1)C=CC=CC2(=CC=CC=C2)'
    mol = Chem.MolFromSmiles(smiles)

    molecule_nx = mol_to_nx(mol)
    #molecule_nx.add_node('a', num = 'a', pos=(1.8096858188806708, -0.2))
    #molecule_nx.add_edge('a', 9)

    molecule_atom = nx.get_node_attributes(molecule_nx, 'num')
    pos = nx.get_node_attributes(molecule_nx, 'pos')

    nx.draw(molecule_nx,
            pos,
            labels=molecule_atom,
            with_labels=True,
            node_color='white',
            connectionstyle='bar, fraction=0.03',
            arrowstyle='-')
    plt.show()

# adjacency matrix
'''matrix = nx.to_numpy_matrix(molecule_nx)
print(matrix)'''
