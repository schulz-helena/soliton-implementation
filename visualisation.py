"""Visualisation of a molecule as a graph
"""
import io
import os

import matplotlib.pyplot as plt
import networkx as nx
from PIL import Image

from soliton_graph import SolitonGraph


class Visualisation:

    @staticmethod
    def visualize_soliton_graph(soliton_graph: SolitonGraph, bindings: dict, show: bool, to_image: bool): # bindings as an extra argument so we can use this method in animation (there we have different bindings for each time step)
        """Plot a visualisation of a soliton graph

        Args:
            soliton_graph (SolitonGraph): soliton graph that should be visualized
            bindings (dict): binding types for each edge in the graph
            show (bool): if the plot should be displayed or not (use False when using visualisation in animation, use True when using it as a stand-alone visualisation)
            title (str): title the saved image of the visualisation should have (if title == None then image is not saved)
        """

        def plot_double_edges(graph: nx.Graph, bindings: dict, double_edge_positions: dict):
            """Plot a second edge for each edge that is a double edhe

            Args:
                graph (nx.Graph): graph the edges are in
                bindings (dict): binding types for each edge in the graph
                double_edge_positions (dict): positions for a second line for each edge in the graph
            """
            for edge in graph.edges:
                if bindings[edge] == 2:
                    x_values, y_values = double_edge_positions[edge]
                    plt.plot(x_values, y_values, color = 'black', linewidth = 1)

        # remove current graph visualisation picture from database folder
        #if os.path.exists('database/graph.jpg'):
            #os.remove('database/graph.jpg')
        plt.clf() # clear plot
        #labels = nx.get_node_attributes(soliton_graph.graph, 'label')
        pos = nx.get_node_attributes(soliton_graph.graph, 'pos')
        if soliton_graph.graph.number_of_nodes() >= 60:
            node_size = 80
            font_size = 5
        elif soliton_graph.graph.number_of_nodes() >= 40:
            node_size = 100
            font_size = 7
        elif (soliton_graph.graph.number_of_nodes() - len(soliton_graph.exterior_nodes)) > 26: # if node labels consist of 2 chars ("aa", ...)
            node_size = 150
            font_size = 9
        else: # use default values
            node_size = 300
            font_size = 12

        plt.axis('equal')
        plot_double_edges(soliton_graph.graph, bindings, soliton_graph.double_edge_positions)
        nx.draw(soliton_graph.graph,
                pos,
                labels=soliton_graph.labels,
                with_labels=True,
                node_color='white',
                node_size = node_size,
                font_size = font_size)
        if to_image == True:        
            #plt.savefig(f'database/{title}.jpg', bbox_inches='tight', format='jpg', dpi=1200)
            # return PIL Image instead of plot
            buf = io.BytesIO()
            plt.savefig(buf, bbox_inches='tight', format='jpg', dpi=1200)
            buf.seek(0)
            im = Image.open(buf)
            im = im.convert("RGBA")
            buf.close()
            return im
        if show == True:
            plt.show()

if __name__ == "__main__":
    
    my_graph = SolitonGraph('C1=CC=CC=C1C{1}=CC{2}=CC=CC2=CC=CC=C2')
    my_graph.validate_soliton_graph()
    Visualisation.visualize_soliton_graph(my_graph, my_graph.bindings, True, None)
    
