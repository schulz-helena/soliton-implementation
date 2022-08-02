"""Animation of a soliton traversing a graph
"""
import io

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import animation
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from PIL import Image

from soliton_graph import SolitonGraph
from soliton_path import SolitonPath
from visualisation import Visualisation


class Animation:

    def __init__(self, soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        self.soliton_graph = soliton_graph
        self.soliton_path = soliton_path
        self.plots_and_arrays = self.list_of_plots_and_arrays()

    def graph_animation(self):
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
            y = [position[1]-0.05]
    
            Visualisation.visualize_soliton_graph(soliton_graph, soliton_path.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black") # plot soliton on top


        fig, ax = plt.subplots() # Build plot
        frames = len(self.soliton_path.path) # as much frames as there are nodes in the soliton path
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(ax, self.soliton_graph, self.soliton_path)) # interval in milliseconds(default 200)

        return ani
        #ani.save('database/animation.gif', writer='ffmpeg')

    def list_of_plots_and_arrays(self):
        plots_and_arrays = []
        for frame_num in range(0, len(self.soliton_path.path)):
            fig, ax = plt.subplots(tight_layout = True, dpi = 500) # no unwanted white spaces and resolution of 200
            canvas = FigureCanvas(fig) # necessary to convert into numpy array later
            ax.clear()
            ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

            pos = nx.get_node_attributes(self.soliton_graph.graph, 'pos')
            node = self.soliton_path.path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05] #TODO: implement smart algorithm that computes best positons and markersize for soliton pebble
    
            Visualisation.visualize_soliton_graph(self.soliton_graph, self.soliton_path.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black") # plot soliton on top

            # convert plot into numpy array
            canvas.draw()
            array_image = np.frombuffer(canvas.tostring_rgb(), dtype="uint8")
            array_image = array_image.reshape(canvas.get_width_height()[::-1] + (3,))

            plots_and_arrays.append((fig, array_image))
            plt.close(fig) # otherwise all created figures are retained in memory -> causes segmentation fault
        return plots_and_arrays

    def list_of_pil_images(self):
        pil_images = []
        for element in self.plots_and_arrays:
            array_image = element[1]
            im = Image.fromarray(np.uint8(array_image))
            pil_images.append(im)
        return pil_images

    '''def graph_animation(self):
        """Animation function that uses update function and saves animation

        Args:
            soliton_graph (SolitonGraph): graph that should be animated
            soliton_path (SolitonPath): path the soliton should traverse in the animation
        """

        def update(frame_num: int, fig: Figure, ax: Axes, plots_and_arrays: list):
            """Update function that plots graph while soliton animation

            Args:
                frame_num (int): frame number (is increased every time animation calls this function)
                ax (Axes): axes of the plot
                soliton_graph (SolitonGraph): graph that should be animated
                soliton_path (SolitonPath): path the soliton should traverse in the animation
            """

            ax = ax
            fig = plots_and_arrays[frame_num][0]
            plt.plot()

        fig, ax = plt.subplots() # Build plot
        frames = len(self.soliton_path.path) # as much frames as there are nodes in the soliton path
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(fig, ax, self.plots_and_arrays)) # interval in milliseconds(default 200)

        ani.save('database/animation.gif', writer='ffmpeg')
        #return ani'''




if __name__ == "__main__":

    my_graph = SolitonGraph('C1{1}=C{3}C1{=2}')
    my_path = SolitonPath([5, 4, 2, 0, 4, 2, 3], my_graph)
    my_animation = Animation(my_graph, my_path)
    #pil_images = my_animation.list_of_pil_images()
    #for i, im in enumerate(pil_images):
        #im.save(f"database/pic{i}.jpg")
    ani = my_animation.graph_animation()
