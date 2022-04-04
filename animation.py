"""Functions to animate a soliton traversing a graph
"""
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import animation
from matplotlib.axes import Axes


def simple_update(frame_num: int, graph: nx.Graph, ax: Axes, path: list):
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
            position = pos[node]
            x = position[0]
            y = position[1]
            ax.add_patch(plt.Circle((x+0.16, y-0.05), radius = 0.04, color='black'))
            if frame_num != 0:
                graph.remove_edge(path[frame_num-1], path[frame_num])
    nx.draw(graph, pos=pos, labels=label, with_labels=True, node_color=color_map, ax=ax)

    # Set the title
    ax.set_title(f"Frame {frame_num}")


def simple_animation(graph: nx.Graph, path: list):
    """Animation function that uses update function and saves animation

    Args:
        graph (nx.Graph): graph that is supposed to be animated
        path (list): path that the soliton should traverse
    """
    # Build plot
    fig, ax = plt.subplots()

    ani = animation.FuncAnimation(fig, simple_update, interval = 400, frames=4, fargs=(graph, ax, path)) #interval in milliseconds(default 200)
    ani.save('test.gif', writer='ffmpeg')

    plt.show()


if __name__ == "__main__":

    g = nx.Graph()
    g.add_node(0, label = '0', pos = (-1.81, 0.02))
    g.add_node(1, label = '1', pos = (-0.72, -0.00))
    g.add_node(2, label = '2', pos = (0.72, 0.00))
    g.add_node(3, label = '3', pos = (1.81, -0.02))
    g.add_edge(0,1)
    g.add_edge(1,2)
    g.add_edge(2,3)

    a_path = [3,2,1,0]

    simple_animation(g, a_path)
