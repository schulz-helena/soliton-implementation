import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import animation
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from PIL import Image
from soliton_classes.soliton_graph import SolitonGraph
from soliton_classes.soliton_path import SolitonPath

from visualisations.visualisation import Visualisation


class Animation:
    """Animation of a soliton (displayed as a black pebble) traversing a graph.
    """

    def __init__(self, soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        """Initializes an `Animation` object by using a soliton graph and a soliton path.
        """
        self.soliton_graph: SolitonGraph = soliton_graph
        """Graph that should be animated."""
        self.soliton_path: SolitonPath = soliton_path
        """Path the soliton should traverse in animation."""
        self.plots_and_arrays: list = self.list_of_plots_and_arrays()
        """Plots and plots as arrays for each timestep in animation."""
        self.pil_images: list = self.list_of_pil_images()
        """`PIL` images for each timestep in animation."""

    def graph_animation(self):
        """Animation method that uses `update` function to build all necessary plots and returns animation.

        Returns:
            animation.FuncAnimation: Animation.
        """

        def update(frame_num: int, ax: Axes, soliton_graph: SolitonGraph, soliton_path: SolitonPath):
            """Update function that plots graph at each timestep.

            Args:
                frame_num (int): Frame number (is increased every time `graph_animation` calls this function).
                ax (Axes): Axes of the plot.
                soliton_graph (SolitonGraph): Graph that should be animated.
                soliton_path (SolitonPath): Path the soliton should traverse in the animation.
            """
            ax.clear()
            ax.axis('equal') # force the x and y axes to have equal number of pixels per data unit (makes circles be round)

            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            node = soliton_path.path[frame_num]
            position = pos[node]
            # position for the black pebble representing the soliton (always located in the bottom right in relation to the node the soliton should be at)
            x = [position[0]+0.12]
            y = [position[1]-0.05]
    
            Visualisation.visualize_soliton_graph(soliton_graph, soliton_path.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black") # plot soliton on top


        fig, ax = plt.subplots() # Build plot
        frames = len(self.soliton_path.path) # as much frames as there are nodes in the soliton path
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(ax, self.soliton_graph, self.soliton_path)) # interval in milliseconds(default 200)

        return ani

    def list_of_plots_and_arrays(self):
        """Builds a list of plots and plots as arrays for each timestep of soliton traversing a path.

        Returns:
            list: `Matplotlib` plots and plots as `Numpy` arrays for each timestep.
        """
        plots_and_arrays = []
        for frame_num, node in enumerate(self.soliton_path.path):
            fig, ax = plt.subplots(tight_layout = True, dpi = 500) # no unwanted white spaces and resolution of 500 (higher resolution can cause segmentation fault)
            canvas = FigureCanvas(fig) # necessary to convert into numpy array later
            ax.clear()
            ax.axis('equal')

            # create plot
            pos = nx.get_node_attributes(self.soliton_graph.graph, 'pos')
            node = self.soliton_path.path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05]
            Visualisation.visualize_soliton_graph(self.soliton_graph, self.soliton_path.bindings_list[frame_num], False, False)
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")

            # convert plot into numpy array
            canvas.draw()
            array_image = np.frombuffer(canvas.tostring_rgb(), dtype="uint8")
            array_image = array_image.reshape(canvas.get_width_height()[::-1] + (3,))

            plots_and_arrays.append((fig, array_image))
            plt.close(fig) # otherwise all created figures are retained in memory -> causes segmentation fault
            
        return plots_and_arrays

    def list_of_pil_images(self):
        """Builds a list of `PIL` images for each timestep of soliton traversing a path.

        Returns:
            list: `PIL` images for each timestep in animation.
        """
        pil_images = []
        for element in self.plots_and_arrays:
            array_image = element[1]
            im = Image.fromarray(np.uint8(array_image))
            pil_images.append(im)

        return pil_images
