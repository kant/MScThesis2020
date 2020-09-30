import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx

from tripsSplice import get_edges_numbers
from tripsSplice import assign_exon_numbers
from tripsSplice import scores_per_exon_number
from tripsSplice import principal_path_in_numbers
from tripsSplice import filter_unsupported_transcripts

''' 
Script to plot any gene as a splice graph in th kamada kawai layout
'''



def nx_plot_spring(gene, sqlite_path_organism, sqlite_path_reads, exclude=True):
    # Conctructs a computational graph using the exons as nodes (splice graph). Exons are labeled with number ID's
    graph = nx.DiGraph()
    supported_transcripts = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)


    edge_numbers = get_edges_numbers(gene, sqlite_path_organism, supported_transcripts, filter=exclude)
    graph.add_edges_from(edge_numbers)
    handles = []
    exon_dictionary = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter=exclude)
    number_of_exons = len(exon_dictionary.keys())

    # Nodes are coloured based on the degree to which they were supported by the reads
    scores = scores_per_exon_number(gene, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=exclude)
    for i in range(number_of_exons + 1):
        if i not in scores:
            scores[i] = 0
    for i in exon_dictionary:
        handles.append(mpatches.Patch(label=str(i) + " : "
                                            + ": " + exon_dictionary[i]))

    connected = list(nx.weakly_connected_components(graph))

    color_lookup = scores
    low, *_, high = sorted(color_lookup.values())
    norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)


    # to improve readability with this layout each connected component is plot individually
    gs = mpl.gridspec.GridSpec(len(connected) + 3, 1)
    plt.gca().set_position(gs[1].get_position(plt.gcf()))
    plt.gca().set_subplotspec(gs[1])


    for i in range(len(connected)):
        subgraph = graph.subgraph(connected[i])
        i_color_lookup = {}
        for j in connected[i]:
            if j not in i_color_lookup:
                i_color_lookup[j] = color_lookup[j]
        subpos = nx.kamada_kawai_layout(subgraph)
        plt.subplot(gs[i])
        nx.draw(subgraph, subpos,
                node_shape='s',
                nodelist=i_color_lookup,
                node_size=150,
                node_color=[mapper.to_rgba(i) for i in i_color_lookup.values()],
                with_labels=True)

    # Include blank subplots to make room for the legend
    blank = nx.Graph()
    plt.subplot(gs[len(connected)])
    nx.draw(blank)
    plt.subplot(gs[len(connected) + 1])
    nx.draw(blank)
    plt.subplot(gs[len(connected) + 2])
    nx.draw(blank)

    # plt.legend(loc="lower left", prop={"size": 4}, handles=handles)
    title = 'Kamada Kawai Layout Network of ' + gene + ' Exons'
    plt.figtext(0.0, 0.95, title)
    principal_path = "Principal Path: ", principal_path_in_numbers(gene, sqlite_path_organism, supported_transcripts, filter=exclude)
    plt.figtext(0.65,0.95,principal_path)
    plt.margins(x=0, y=2)
    plt.show()

gene = "phpt1"
transcript = "ENST00000463215"
sqlite_path_organism = "homo_sapiens.v2.sqlite"
sqlite_path_reads = "SRR2433794.sqlite"


nx_plot_spring("phpt1", sqlite_path_organism, sqlite_path_reads, exclude=False)