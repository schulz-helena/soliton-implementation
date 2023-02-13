"""Animation of a soliton traversing a graph.
"""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import animation
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.colors import hsv_to_rgb
from PIL import Image

from soliton_automata.soliton_classes.soliton_graph import SolitonGraph
from soliton_automata.soliton_classes.soliton_path import SolitonPath
from soliton_automata.soliton_classes.traversal import Traversal
from soliton_automata.visualisations.visualisation import Visualisation


class Animation:
    """Animations of a soliton/ solitons (displayed as pebbles) traversing a graph.
    """

    @staticmethod
    def graph_animation(soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        """Animation of a single soliton traversing a graph. 
        Uses `update` function to build all necessary plots and returns animation.

        Args:
            soliton_graph (SolitonGraph): The graph that is traversed.
            soliton_path (SolitonPath): Path the soliton should traverse in the animation.
        Returns:
            animation.FuncAnimation: Animation.
        """

        def update(frame_num: int, ax: Axes, soliton_graph: SolitonGraph, soliton_path: SolitonPath):
            """Update function that plots graph at each timestep.

            Args:
                frame_num (int): Frame number (is increased every time `graph_animation` calls this function).
                ax (Axes): Axes of the plot.
                soliton_graph (SolitonGraph): The graph that is traversed.
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
        frames = len(soliton_path.path) # as much frames as there are nodes in the soliton path
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(ax, soliton_graph, soliton_path)) # interval in milliseconds(default 200)

        return ani


    @staticmethod
    def graph_animation_multiwave(soliton_graph: SolitonGraph, traversal: Traversal):
        """Animation of multiple solitons traversing a graph. 
        Uses `update` function to build all necessary plots and returns animation.

        Args:
            soliton_graph (SolitonGraph): The graph that is traversed.
            traversal (Traversal): The traversal that contains bindings and soliton positions for each timestep.
        Returns:
            animation.FuncAnimation: Animation.
        """

        def update(frame_num: int, ax: Axes, soliton_graph: SolitonGraph, traversal: Traversal, color_margin: int):
            """Update function that plots graph at each timestep.

            Args:
                frame_num (int): Frame number (is increased every time `graph_animation` calls this function).
                ax (Axes): Axes of the plot.
                soliton_graph (SolitonGraph): The graph that is traversed.
                traversal (Traversal): The traversal that contains bindings and soliton positions for each timestep.
                color_margin(int): Margin between the hsv colors of the different solitons.
            """
            ax.clear()
            ax.axis('equal')

            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            soliton_positions = {}
            unique_positions = []
            positions = traversal.pos[frame_num]
            for soliton in positions:
                node = positions[soliton]
                if node == -1 or node == -2: # at this timestep the soliton is not in the graph
                    soliton_positions[soliton] = 0
                else:
                    position = pos[node]
                    if (position[0]+0.12, position[1]-0.05) not in unique_positions:
                        # put in bottom right corner of node
                        x = position[0]+0.12
                        y = position[1]-0.05
                    else: # if there already is a soliton at this node, put soliton pebble in bottom left corner
                        x = position[0]-0.12
                        y = position[1]-0.05
                    unique_positions.append((x,y))
                    soliton_positions[soliton] = (x,y)
    
            Visualisation.visualize_soliton_graph(soliton_graph, traversal.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            akt_color = 0
            for soliton in soliton_positions:
                if soliton != 1:
                    akt_color += color_margin
                if soliton_positions[soliton] != 0: # if soliton is in the graph at this timestep
                    x = soliton_positions[soliton][0]
                    y = soliton_positions[soliton][1]
                    if soliton == 1:
                        plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")
                    else:
                        rgb = hsv_to_rgb((akt_color/360, 1, 1.0))
                        plt.plot(x, y, marker="o", markersize=6, markeredgecolor=rgb, markerfacecolor=rgb, zorder = 3)


        fig, ax = plt.subplots()
        frames = len(traversal.pos)
        if traversal.soliton_num != 1:
            color_margin = 360 / (traversal.soliton_num - 1)
        else: color_margin = 0
        ani = animation.FuncAnimation(fig, update, interval = 800, frames=frames, fargs=(ax, soliton_graph, traversal, color_margin))

        return ani


    @staticmethod
    def list_of_plots_and_arrays(soliton_graph: SolitonGraph, soliton_path: SolitonPath):
        """Builds a list of plots and plots as arrays for each timestep of soliton traversing a path.

        Args:
            soliton_graph (SolitonGraph): The graph that is traversed.
            soliton_path (SolitonPath): Path the soliton should traverse in the animation.
        Returns:
            list: `Matplotlib` plots and plots as `Numpy` arrays for each timestep.
        """
        plots_and_arrays = []
        for frame_num, node in enumerate(soliton_path.path):
            fig, ax = plt.subplots(tight_layout = True, dpi = 500) # no unwanted white spaces and resolution of 500 (higher resolution can cause segmentation fault)
            canvas = FigureCanvas(fig) # necessary to convert into numpy array later
            ax.clear()
            ax.axis('equal')

            # create plot
            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            node = soliton_path.path[frame_num]
            position = pos[node]
            x = [position[0]+0.12]
            y = [position[1]-0.05]
            Visualisation.visualize_soliton_graph(soliton_graph, soliton_path.bindings_list[frame_num], False, False)
            plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")

            # convert plot into numpy array
            canvas.draw()
            array_image = np.frombuffer(canvas.tostring_rgb(), dtype="uint8")
            array_image = array_image.reshape(canvas.get_width_height()[::-1] + (3,))

            plots_and_arrays.append((fig, array_image))
            plt.close(fig) # otherwise all created figures are retained in memory -> causes segmentation fault
            
        return plots_and_arrays


    @staticmethod
    def list_of_plots_and_arrays_multiwave(soliton_graph: SolitonGraph, traversal: Traversal):
        """Builds a list of plots and plots as arrays for each timestep of solitons traversing the graph.

        Args:
            soliton_graph (SolitonGraph): The graph that is traversed.
            traversal (Traversal): The traversal that contains bindings and soliton positions for each timestep.
        Returns:
            list: `Matplotlib` plots and plots as `Numpy` arrays for each timestep.
        """
        plots_and_arrays = []
        if traversal.soliton_num != 1:
            color_margin = 360 / (traversal.soliton_num - 1)
        else: color_margin = 0
        for frame_num, node in enumerate(traversal.pos):
            fig, ax = plt.subplots(tight_layout = True, dpi = 500) # no unwanted white spaces and resolution of 500 (higher resolution can cause segmentation fault)
            canvas = FigureCanvas(fig) # necessary to convert into numpy array later
            ax.clear()
            ax.axis('equal')

            # create plot
            pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
            soliton_positions = {}
            unique_positions = []
            positions = traversal.pos[frame_num]
            for soliton in positions:
                node = positions[soliton]
                if node == -1 or node == -2: # at this timestep the soliton is not in the graph
                    soliton_positions[soliton] = 0
                else:
                    position = pos[node]
                    if (position[0]+0.12, position[1]-0.05) not in unique_positions:
                        # put in bottom right corner of node
                        x = position[0]+0.12
                        y = position[1]-0.05
                    else: # if there already is a soliton at this node, put soliton pebble in bottom left corner
                        x = position[0]-0.12
                        y = position[1]-0.05
                    unique_positions.append((x,y))
                    soliton_positions[soliton] = (x,y)
    
            Visualisation.visualize_soliton_graph(soliton_graph, traversal.bindings_list[frame_num], False, False) # use visualisation of graph at current timestep
            akt_color = 0
            for soliton in soliton_positions:
                if soliton != 1:
                    akt_color += color_margin
                if soliton_positions[soliton] != 0: # if soliton is in the graph at this timestep
                    x = soliton_positions[soliton][0]
                    y = soliton_positions[soliton][1]
                    if soliton == 1:
                        plt.plot(x, y, marker="o", markersize=6, markeredgecolor="black", markerfacecolor="black")
                    else:
                        rgb = hsv_to_rgb((akt_color/360, 1, 1.0))
                        plt.plot(x, y, marker="o", markersize=6, markeredgecolor=rgb, markerfacecolor=rgb, zorder = 3)

            # convert plot into numpy array
            canvas.draw()
            array_image = np.frombuffer(canvas.tostring_rgb(), dtype="uint8")
            array_image = array_image.reshape(canvas.get_width_height()[::-1] + (3,))

            plots_and_arrays.append((fig, array_image))
            plt.close(fig) # otherwise all created figures are retained in memory -> causes segmentation fault
            
        return plots_and_arrays


    @staticmethod
    def list_of_pil_images(plots_and_arrays: list):
        """Builds a list of `PIL` images for each timestep of soliton traversing a path/ solitons traversing the graph.

        Args:
            plots_and_arrays (list): `Matplotlib` plots and plots as `Numpy` arrays for each timestep. 
        Returns:
            list: `PIL` images for each timestep in animation.
        """
        pil_images = []
        for element in plots_and_arrays:
            array_image = element[1]
            im = Image.fromarray(np.uint8(array_image))
            pil_images.append(im)

        return pil_images
