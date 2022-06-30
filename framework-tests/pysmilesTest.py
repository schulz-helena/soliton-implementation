import matplotlib.pyplot as plt
import networkx as nx
from pysmiles import read_smiles

#smiles = 'C1CC[13CH2]CC1C1CCCCC1'
smiles = 'C1=CC=CC=C1C=CC=CC2=CC=CC=C2'
mol = read_smiles(smiles, reinterpret_aromatic=False)
#print(mol)

#print(mol.nodes(data='element'))
order = nx.get_edge_attributes(mol, 'order')
print(order)
if (12,13) in order:
    print(order[(12,13)])
#nx.draw(mol)
#plt.show()
