"""Functions to animate a soliton traversing a graph
"""
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import animation
from matplotlib.axes import Axes
from pysmiles import read_smiles

import visualisation as v


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
    #ax.set_title(f"Frame {frame_num}") #set the title
    for node in graph:
        if node == path[frame_num]:
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
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05] #TODO: implement smart algorithm that computes best positons and markersize for soliton pebble
    
    for edge in graph.edges:
        if bindings[edge] == 2:
            x_values, y_values = xs_and_ys[edge]
            plt.plot(x_values, y_values, color = 'black', linewidth = 1)
    
    nx.draw(graph, pos=pos, labels=label, with_labels=True, node_color = 'white', ax=ax)
    plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")



def simple_animation(graph: nx.Graph, path: list, bindings: dict, xs_and_ys: dict, frames: int):
    """Animation function that uses update function and saves animation

    Args:
        graph (nx.Graph): graph that is supposed to be animated
        path (list): path that the soliton should traverse
    """
    # Build plot
    fig, ax = plt.subplots()

    ani = animation.FuncAnimation(fig, simple_update, interval = 800, frames=frames, fargs=(graph, ax, path, bindings, xs_and_ys)) #interval in milliseconds(default 200)
    ani.save('test.gif', writer='ffmpeg')

    #plt.show()
    return bindings


if __name__ == "__main__":

    #smiles = 'C1{1}=C{2}C1{=3}'
    #a_path = [1,0,2,4,5]
    smiles = 'C1=CC{1}=CC=C1C{=2}=CC=C{=3}C2=C{4}C=CC=C2'
    a_path = [3,2,4,5,6,7,8]

    ext_nodes, smi, rdkit_smi = v.transform_user_input(smiles)
    mol_pysmiles = read_smiles(smi, reinterpret_aromatic=False)
    binds = nx.get_edge_attributes(mol_pysmiles, 'order')
    print(binds)

    g, nix = v.mol_to_nx3(ext_nodes, smi, rdkit_smi)
    xs_and_ys = v.double_edge_positions_dict(g)


    bind = simple_animation(g, a_path, binds, xs_and_ys, len(a_path))

    print(bind)
