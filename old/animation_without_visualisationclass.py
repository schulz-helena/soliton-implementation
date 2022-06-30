"""Functions to animate a soliton traversing a graph
"""
from dataclasses import dataclass

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import animation
from matplotlib.axes import Axes

from soliton_graph import SolitonGraph
from soliton_path import SolitonPath


@dataclass
class Animation:

    soliton_graph: SolitonGraph
    soliton_path: SolitonPath

    def simple_animation(self, frames: int):
        """Animation function that uses update function and saves animation

        Args:
            frames (int): number of frames for animation
        """

        def simple_update(frame_num: int, ax: Axes, graph: nx.Graph, path: list, bindings_list: list, double_edge_positions: dict):
            """Update function that plots graph while soliton animation

            Args:
                frame_num (int): frame number (is increased every time animation calls this function)
                ax (Axes): axes of the plot
            """
            ax.clear()
            ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

            label = nx.get_node_attributes(graph, 'label')
            pos = nx.get_node_attributes(graph, 'pos')
            #ax.set_title(f"Frame {frame_num}") #set the title
            node = path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05] #TODO: implement smart algorithm that computes best positons and markersize for soliton pebble
    
            for edge in bindings_list[frame_num]:
                if bindings_list[frame_num][edge] == 2:
                    x_values, y_values = double_edge_positions[edge]
                    plt.plot(x_values, y_values, color = 'black', linewidth = 1)
    
            nx.draw(graph, pos=pos, labels=label, with_labels=True, node_color = 'white', ax=ax)
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")


        fig, ax = plt.subplots() # Build plot
        ani = animation.FuncAnimation(fig, simple_update, interval = 800, frames=frames, fargs=(ax, self.soliton_graph.graph, self.soliton_path.path, self.soliton_path.bindings_list, self.soliton_graph.double_edge_positions)) #interval in milliseconds(default 200)
        ani.save('test.gif', writer='ffmpeg')
        #plt.show()

if __name__ == "__main__":

    my_graph = SolitonGraph('C1=CC=CC=C1C=C{1}C=CC{2}=CC2=CC=CC=C2')
    my_path = SolitonPath([8,7,6,5,4,3,2,1,0,5,6,7,9,10,11,12], my_graph)
    my_animation = Animation(my_graph, my_path)

    my_animation.simple_animation(len(my_path.path))
