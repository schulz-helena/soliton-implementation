"""Visualisation of a molecule as a graph
"""
import matplotlib.pyplot as plt
import networkx as nx

from soliton_graph import SolitonGraph


class Visualisation:

    @staticmethod
    def visualize_soliton_graph(soliton_graph: SolitonGraph, bindings: dict, show: bool): # bindings as an extra argument so we can use this method in animation (there we have different bindings for each time step)
        """Plot a visualisation of a soliton graph

        Args:
            soliton_graph (SolitonGraph): soliton graph that should be visualized
            bindings (dict): binding types for each edge in the graph
            show (bool): if the plot should be displayed or not (use False when using visualisation in animation, use True when using it as a stand-alone visualisation)
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

        labels = nx.get_node_attributes(soliton_graph.graph, 'label')
        pos = nx.get_node_attributes(soliton_graph.graph, 'pos')

        plt.axis('equal')
        plot_double_edges(soliton_graph.graph, bindings, soliton_graph.double_edge_positions)
        nx.draw(soliton_graph.graph,
                pos,
                labels=labels,
                with_labels=True,
                node_color='white')
        if show == True:
            plt.show()

if __name__ == "__main__":
    
    my_graph = SolitonGraph('C1=CC=CC=C1C{1}=CC{2}=CC=CC2=CC=CC=C2')
    my_graph.validate_soliton_graph()
    Visualisation.visualize_soliton_graph(my_graph, my_graph.bindings, True)
