"""Animation of a soliton traversing a graph
"""
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import animation
from matplotlib.axes import Axes

from soliton_graph import SolitonGraph
from soliton_path import SolitonPath
from visualisation import Visualisation


class Animation:

    @staticmethod
    def graph_animation(soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        """Animation function that uses update function and saves animation

        Args:
            frames (int): number of frames for animation
        """

        def update(frame_num: int, ax: Axes, soliton_graph: SolitonGraph, soliton_path: SolitonPath):
            """Update function that plots graph while soliton animation

            Args:
                frame_num (int): frame number (is increased every time animation calls this function)
                ax (Axes): axes of the plot
            """
            ax.clear()
            ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            node = soliton_path.path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05] #TODO: implement smart algorithm that computes best positons and markersize for soliton pebble
    
            Visualisation.visualize_soliton_graph(soliton_graph, soliton_path.bindings_list[frame_num], False)
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")


        fig, ax = plt.subplots() # Build plot
        frames = len(soliton_path.path) #as much frames as there are nodes in the soliton path
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(ax, soliton_graph, soliton_path)) #interval in milliseconds(default 200)
        ani.save('animation.gif', writer='ffmpeg')

if __name__ == "__main__":

    my_graph = SolitonGraph('C1=CC=CC=C1C=C{1}C=CC{2}=CC2=CC=CC=C2')
    my_path = SolitonPath([8,7,6,5,4,3,2,1,0,5,6,7,9,10,11,12], my_graph)
    Animation.graph_animation(my_graph, my_path)
