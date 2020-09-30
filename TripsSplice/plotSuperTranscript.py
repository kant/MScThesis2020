import matplotlib.pyplot as plt

from tripsSplice import get_transcript_modules
from tripsSplice import get_reads_per_genomic_location_asite
from tripsSplice import get_module_pileup
from tripsSplice import supertranscript_positions
from tripsSplice import supertranscript_colors
from tripsSplice import supertranscript_widths_starts
from tripsSplice import supertranscript_edgecolors
from tripsSplice import get_principal_modules
from tripsSplice import filter_unsupported_transcripts

gene = "phpt1"
transcript = "ENST00000463215"
sqlite_path_organism = "homo_sapiens.v2.sqlite"
sqlite_path_reads = ["SRR2433794.sqlite"]

def plot_supertranscript(gene, sqlite_path_organism, sqlite_path_reads, exclude):

    supported_transcripts = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)

    modules = get_transcript_modules(gene, sqlite_path_organism, sqlite_path_reads, filter=exclude)
    genomic_read_dictionary = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                             supported_transcripts, filter=exclude)
    scores = get_module_pileup(modules, genomic_read_dictionary)
    fig, ax = plt.subplots()
    positions = supertranscript_positions(modules)
    ax.broken_barh(positions, (0, 1), facecolors = 'blue', edgecolors = 'black')
    ax.set_ylim(0, 15)
    ax.set_xlim(0, positions[-1][0] + positions[-1][1])
    plt.show()

# plot_supertranscript(gene, sqlite_path_organism, sqlite_path_reads)

def plot_supertranscript_barh(gene, sqlite_path_organism, sqlite_path_reads, exclude=True):
    # horizontal bar chart method of plotting a supertranscript. A modules support is
    # reflected in the color concentration.
    supported_transcripts = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)

    try:
        modules = get_transcript_modules(gene, sqlite_path_organism, supported_transcripts, filter=exclude)
    except ValueError:
        print("There are no transcripts where all exons have supporting reads for this gene.")
        print()
        user_input = input("Set exclude to False? (True/False): ")
        if (user_input == "True") or (user_input == "T"):
            exclude = False
            modules = get_transcript_modules(gene, sqlite_path_organism, supported_transcripts, filter=False)
        else:
            print()
            print("No fully supported transcripts. No unfiltered plot requested")
            return

    genomic_read_dictionary = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                             supported_transcripts, filter=exclude)
    scores = get_module_pileup(modules, genomic_read_dictionary)
    colors = supertranscript_colors(scores)
    widths, starts = supertranscript_widths_starts(modules)
    fig, ax = plt.subplots(figsize=(10,1.5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(False)
    c = plt.get_cmap('Blues')(colors)

    principal_modules = get_principal_modules(gene, sqlite_path_organism, modules, supported_transcripts, exclude)
    edgecolors = supertranscript_edgecolors(modules, principal_modules)

    ax.barh(1, widths, left=starts, height=0.5, edgecolor=edgecolors, color=c)
    plt.show()

plot_supertranscript_barh('brca1', sqlite_path_organism, sqlite_path_reads, exclude=False)

