"""
script to plot a splice graph of a gene with specified layout based on genomic
position of exons.
"""
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx

from tripsSplice import exon_positions_with_numbers
from tripsSplice import get_y_values
from tripsSplice import get_edges_numbers
from tripsSplice import assign_exon_numbers
from tripsSplice import scores_per_exon_number
from tripsSplice import principal_path_in_numbers
from tripsSplice import filter_unsupported_transcripts
from tripsSplice import get_edges_scores_numbers
from tripsSplice import get_scores_for_gene
from tripsSplice import relabel_nodes_L_R

from timer import Timer

def plot_genomic_pos(gene, sqlite_path_organism, sqlite_path_reads, exclude=True):
    # Conctructs a computational graph using the exons as nodes (splice graph). Exons are labeled with number ID's
    # exon_positions_with_numbers gets the genomic position of the exon start for each node.
    # get y_values sacks nodes up the y axis depending how populated that locus is with the annotated exons.
    # this ensures no overlapping of nodes
    t = Timer()
    t.start()

    supported_transcripts = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)
    try:
        numbers_position = exon_positions_with_numbers(gene, sqlite_path_organism, supported_transcripts, exclude)
    except ValueError:
        print("There are no transcripts where all exons have supporting reads for this gene.")
        print()
        user_input = input("Set exclude to False? (True/False): ")
        if (user_input == "True") or (user_input == "T"):
            exclude = False
            numbers_position = exon_positions_with_numbers(gene, sqlite_path_organism, supported_transcripts, False)
        else:
            print()
            print("No fully supported transcripts. No unfiltered plot requested")
            return

    number_x_y = get_y_values(numbers_position)
    node_labels = relabel_nodes_L_R(number_x_y)
    graph = nx.DiGraph()
    for node in number_x_y:
        graph.add_node(node, pos=(number_x_y[node][0], number_x_y[node][1]))

    edge_numbers = get_edges_numbers(gene, sqlite_path_organism, supported_transcripts, exclude)
    graph.add_edges_from(edge_numbers, color="black")

    handles = []
    exon_dict = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, exclude)
    number_of_exons = len(exon_dict.keys())

    # Nodes are coloured based on the degree to which they were supported by the reads
    scores, junction_pileup = get_scores_for_gene(gene, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=True)
    for i in range(1, number_of_exons + 1):
        if i not in scores:
            scores[i] = 0
    principal_path = principal_path_in_numbers(gene, sqlite_path_organism, supported_transcripts, exclude)
    if exclude:
        zero_support = []
        for exon in scores:
            if (scores[exon] == 0) and (exon not in principal_path):
                graph.remove_node(exon)
                zero_support.append(exon)

        for exon in zero_support:
            del scores[exon]

    for i in exon_dict:
        handles.append(mpatches.Patch(label=str(i) + " : " + ": " + exon_dict[i]))

    color_lookup = scores
    low, *_, high = sorted(color_lookup.values())
    norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)

    pos = nx.get_node_attributes(graph, 'pos')

    # Path of principal isoform as per APPRIS database is coloured red
    principal_edges = []
    edges = graph.edges()

    for i in range(len(principal_path) - 1):
        if principal_path[i] not in list(edges):
            graph.add_node(principal_path[i])
        elif principal_path[i + 1] not in list(edges):
            graph.add_node(principal_path[i + 1])

        principal_edges.append((principal_path[i], principal_path[i + 1]))

    graph.add_edges_from(principal_edges, color="r")
    colors = [graph[u][v]['color'] for u, v in edges]

    edge_labels = get_edges_scores_numbers(gene, sqlite_path_organism, supported_transcripts, junction_pileup, filter=exclude)
    nx.draw(graph, pos,
            edges=edges,
            edge_color=colors,
            node_shape='s',
            nodelist=color_lookup,
            node_size=150,
            node_color=[mapper.to_rgba(i) for i in color_lookup.values()],
            with_labels=False,
            linewidths=1,
            font_size=8)

    ax = plt.gca()  # to get the current axis
    ax.collections[0].set_edgecolor("#000000")

    # nx.draw_networkx_labels(graph, pos, node_labels, font_size=7) # Node labels from L to R
    nx.draw_networkx_labels(graph, pos, scores, font_size=5) #Node lables as scores per that exon
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_color = 'blue', font_size=6)
    principal_path_new = [node_labels[i] for i in principal_path]
    principal_text = "Principal Path: ", principal_path_new

    # plt.figtext(0.01, 0.95, principal_text)
    plt.margins(x=0, y=0.5)
    t.stop()
    plt.show()


gene = "phpt1"
transcript = "ENST00000463215"
sqlite_path_organism = "homo_sapiens.v2.sqlite"
sqlite_path_reads = ["SRR2433794.sqlite"]

# plot_genomic_pos("brca1", sqlite_path_organism, sqlite_path_reads, exclude=True)

plot_genomic_pos("phpt1", sqlite_path_organism, sqlite_path_reads, exclude=True)


# plot_genomic_pos(gene, sqlite_path_organism, sqlite_path_reads, exclude=False)

# plot_genomic_pos('igf2', sqlite_path_organism, sqlite_path_reads, exclude=False)
