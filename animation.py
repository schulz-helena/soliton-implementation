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
            #position = pos[node]
            #x = position[0]
            #y = position[1]
            #ax.add_patch(plt.Circle((x+0.16, y-0.05), radius = 0.04, color='black'))
            #if frame_num != 0:
                #graph.remove_edge(path[frame_num-1], path[frame_num])
    
    for edge in graph.edges:
        if bindings[edge] == 2:
            x_values, y_values = xs_and_ys[edge]
            plt.plot(x_values, y_values, color = 'black', linewidth = 1)

    nx.draw(graph, pos=pos, labels=label, with_labels=True, node_color = color_map, ax=ax)

    # Set the title
    ax.set_title(f"Frame {frame_num}")


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

    '''g = nx.Graph()
    g.add_node(0, label = '0', pos = (-1.81, 0.02))
    g.add_node(1, label = '1', pos = (-0.72, -0.00))
    g.add_node(2, label = '2', pos = (0.72, 0.00))
    g.add_node(3, label = '3', pos = (1.81, -0.02))
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(2,3)'''

    smiles = 'C1{1}=C{3}C1{=2}'

    ext_nodes, smi, rdkit_smi = v.transform_user_input(smiles)
    mol_pysmiles = read_smiles(smi, reinterpret_aromatic=False)
    binds = nx.get_edge_attributes(mol_pysmiles, 'order')
    print(binds)

    g, nix = v.mol_to_nx3(ext_nodes, smi, rdkit_smi)
    xs_and_ys = v.double_edge_positions_dict(g)

    a_path = [1,0,2,4,5]

    bind = simple_animation(g, a_path, binds, xs_and_ys, len(a_path))

    print(bind)
