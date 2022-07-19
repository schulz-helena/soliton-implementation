"""Animation of a soliton traversing a graph
"""
import io

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import animation
from matplotlib.axes import Axes
from PIL import Image

from soliton_graph import SolitonGraph
from soliton_path import SolitonPath
from visualisation import Visualisation


class Animation:

    @staticmethod
    def graph_animation(soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        """Animation function that uses update function and saves animation

        Args:
            soliton_graph (SolitonGraph): graph that should be animated
            soliton_path (SolitonPath): path the soliton should traverse in the animation
        """

        def update(frame_num: int, ax: Axes, soliton_graph: SolitonGraph, soliton_path: SolitonPath):
            """Update function that plots graph while soliton animation

            Args:
                frame_num (int): frame number (is increased every time animation calls this function)
                ax (Axes): axes of the plot
                soliton_graph (SolitonGraph): graph that should be animated
                soliton_path (SolitonPath): path the soliton should traverse in the animation
            """
            ax.clear()
            ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            node = soliton_path.path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05] #TODO: implement smart algorithm that computes best positons and markersize for soliton pebble
    
            Visualisation.visualize_soliton_graph(soliton_graph, soliton_path.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black") # plot soliton on top


        fig, ax = plt.subplots() # Build plot
        frames = len(soliton_path.path) # as much frames as there are nodes in the soliton path
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(ax, soliton_graph, soliton_path)) # interval in milliseconds(default 200)

        #ani.save('database/animation.gif', writer='ffmpeg')
        return ani

    @staticmethod
    def list_of_plots(soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        list_of_plots = []
        for frame_num in range(0, len(soliton_path.path)):
            fig, ax = plt.subplots()
            ax.clear()
            ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            node = soliton_path.path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05] #TODO: implement smart algorithm that computes best positons and markersize for soliton pebble
    
            Visualisation.visualize_soliton_graph(soliton_graph, soliton_path.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black") # plot soliton on top
            list_of_plots.append(plt)
        return list_of_plots

    @staticmethod
    def list_of_pil_images(list_of_plots: list):
        list_of_pil_images = []
        for plt in list_of_plots:
            buf = io.BytesIO()
            plt.savefig(buf, bbox_inches='tight', format='jpg', dpi=1200)
            buf.seek(0)
            im = Image.open(buf)
            im = im.convert("RGBA")
            list_of_pil_images.append(im)
            buf.flush()
            buf.close()
        return list_of_pil_images




if __name__ == "__main__":

    my_graph = SolitonGraph('C1{1}=C{3}C1{=2}')
    my_path = SolitonPath([5, 4, 2, 0, 4, 2, 3], my_graph)
    plots = Animation.list_of_plots(my_graph, my_path)
    images = Animation.list_of_pil_images(plots)
    for im in images:
        im.show()
    #ani = Animation.graph_animation(my_graph, my_path)
    #ani.save('database/animation.gif', writer='ffmpeg')
