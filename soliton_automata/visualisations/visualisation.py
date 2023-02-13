"""Visualisation of a soliton graph.
"""
import io

import matplotlib.pyplot as plt
import networkx as nx
from PIL import Image
from mini_soliton_automata.soliton_classes.soliton_graph import SolitonGraph


class Visualisation:
    """Visualisation of a molecule as a graph.
    """

    @staticmethod
    def visualize_soliton_graph(soliton_graph: SolitonGraph, bindings: dict, show: bool, to_image: bool): # bindings as an extra argument so we can use this method in animation (there we have different bindings for each time step)
        """Plots a visualisation of a soliton graph.

        Args:
            soliton_graph (SolitonGraph): Soliton graph that should be visualised.
            bindings (dict): Binding types for each edge in the graph.
            show (bool): Whether or not the plot should be displayed (use `False` when using visualisation in animation, use `True` when using it as a stand-alone visualisation).
            to_image (bool): Whether or not the plot should be returned as a `PIL` image.
        """

        def plot_double_edges(graph: nx.Graph, bindings: dict, double_edge_positions: dict):
            """Plots a second edge for each edge that is a double edge.

            Args:
                graph (nx.Graph): Graph the edges are contained in.
                bindings (dict): Binding types for each edge in the graph.
                double_edge_positions (dict): Positions for a second line for each edge in the graph.
            """
            for edge in graph.edges:
                if bindings[edge] == 2:
                    x_values, y_values = double_edge_positions[edge]
                    plt.plot(x_values, y_values, color = 'black', linewidth = 1)

        plt.clf() # clear plot
        pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
        # use different node and font sizes for different numbers of nodes
        if soliton_graph.graph.number_of_nodes() >= 60:
            node_size = 80
            font_size = 5
        elif soliton_graph.graph.number_of_nodes() >= 40:
            node_size = 100
            font_size = 7
        elif (soliton_graph.graph.number_of_nodes() - len(soliton_graph.exterior_nodes)) > 26: # if node labels consist of 2 chars ("aa", ...)
            node_size = 150
            font_size = 9
        else: # use (networkx's) default values
            node_size = 300
            font_size = 12
        plt.axis('equal')
        # plot double edges first and then plot rest of the graph on top of it 
        plot_double_edges(soliton_graph.graph, bindings, soliton_graph.double_edge_positions)
        nx.draw(soliton_graph.graph,
                pos,
                labels=soliton_graph.labels,
                with_labels=True,
                node_color='white',
                node_size = node_size,
                font_size = font_size)
                
        if to_image is True:
            # return PIL Image instead of plot
            buf = io.BytesIO()
            plt.savefig(buf, bbox_inches='tight', format='jpg', dpi=1200)
            buf.seek(0)
            im = Image.open(buf)
            im = im.convert("RGBA")
            buf.close()
            return im
        if show is True:
            plt.show()
