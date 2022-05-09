"""Visualisation of a molecule as a graph
"""
import matplotlib.pyplot as plt
import networkx as nx

from soliton_graph import SolitonGraph


class Visualisation:

    @staticmethod
    def visualize_soliton_graph(soliton_graph: SolitonGraph, bindings: dict, show: bool): #bindings as an extra argument so we can use this method in animation (there we have different bidnings for each time step)
        """Plot a visualisation of a soliton graph

        Args:
            soliton_graph (SolitonGraph): soliton graph that should be visualized
        """

        def plot_double_edges(graph: nx.Graph, bindings: dict, double_edge_positions: dict):
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
    my_graph = SolitonGraph('CC=CC')
    my_graph.validate_soliton_graph()
    Visualisation.visualize_soliton_graph(my_graph, my_graph.bindings, True)
