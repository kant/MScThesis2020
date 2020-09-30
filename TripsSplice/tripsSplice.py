import sqlite3

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
from sqlitedict import SqliteDict

gene = "phpt1"
transcript = "ENST00000463215"
sqlite_path_organism = "homo_sapiens.v2.sqlite"
sqlite_path_reads = ["SRR2433794.sqlite"]


def get_gene_info(gene_id, sqlite_path_organism):
    # get the transcript ID, exon junctions and transcript sequence from the specified sqlite file for a given gene ID
    gene_id = gene_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript, exon_junctions, sequence FROM transcripts WHERE gene = (:gene);', {'gene': gene_id})
    gene_info = c.fetchall()
    return gene_info


def get_genes_principal(gene_id, sqlite_path_organism):
    # get the transcript ID for the given genes principal isoform as per APPRIS
    gene_id = gene_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript FROM transcripts '
              'WHERE gene = (:gene) AND principal = 1;', {'gene': gene_id})
    principal_transcript = c.fetchall()
    return principal_transcript[0][0]


def string_to_integer_list(lst):
    # Return the list of integers from a list of strings. If the list is empty then return 0
    if lst[0] == "":
        lst[0] = "0"
    new_lst = []
    for i in range(len(lst)):
        new_lst.append(int(lst[i]))
    return new_lst


def get_3prime_exon(junction_list, sequence):
    # Part of the process of producing the sequences of all exons in a transcript.
    # this function slices the 3' exon sequence from the transcript sequence returning
    # both the exon sequence and the remaining sequence of the transcript
    cut_off = junction_list[-1]
    exon = sequence[cut_off-1:]
    seq_less_exon = sequence[:cut_off-1]
    return exon, seq_less_exon


def get_transcript_info(transcript_id, sqlite_path_organism):
    # Returns the exon junctions and sequences for the given transcript
    transcript_id = transcript_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT exon_junctions, sequence FROM transcripts WHERE transcript = (:trans);', {'trans': transcript_id})
    transcript_info = c.fetchall()
    return transcript_info


def get_max_min_exon_genomic_positions(exon_info):
    # Of all exons starts return the lowest position and highest position of all stops
    starts = []
    for i, value in enumerate(exon_info):
        starts.append(exon_info[i][1])

    stops = []
    for i, value in enumerate(exon_info):
        stops.append(exon_info[i][2])
    return min(starts), max(stops)


def exons_of_transcript(transcript_id, sqlite_path_organism):
    # For a given transcript return the exon sequences in a list in 5' to 3' direction 
    exon_lst = []
    trans_info = get_transcript_info(transcript_id, sqlite_path_organism)
    exon_junctions = trans_info[0][0].split(",")
    sequence = trans_info[0][1]

    exon_junct_int = string_to_integer_list(exon_junctions)

    while len(exon_junct_int) != 0:
        exon, sequence = get_3prime_exon(exon_junct_int, sequence)
        exon_lst.append(exon)

        exon_junct_int.pop(-1)
    exon_lst.append(sequence)
    return exon_lst[::-1]


def get_exon_coordinate_ranges(sequence, exons, junctions):
    # Return list of transcript coordinate ranges (start position, stop position).
    end = len(sequence)-1
    junctions.append(end)

    ranges = []
    for i in range(len(exons)):
        ranges.append((sequence.find(exons[i]), junctions[i]-1))
    return ranges


def get_exon_info(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # Return the exon starts,stops and transcript for a given gene
    gene = gene.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript, exon_start, exon_stop FROM exons WHERE transcript IN '
              '(SELECT transcript FROM transcripts WHERE gene = (:gene));', {'gene': gene})
    exon_info = c.fetchall()
    if filter:
        supported_exon_info = []
        for exon in exon_info:
            if exon[0] in supported_transcripts:
                supported_exon_info.append(exon)
        return supported_exon_info

    return exon_info


def genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # create a dictionary with transcript_ids as keys and exon coordinates in tuple (start, stop) as values
    # subtract the minumum start codon position of any exon for the gene

    exon_info = get_exon_info(gene, sqlite_path_organism, supported_transcripts, filter=filter)
    minimum, _ = get_max_min_exon_genomic_positions(exon_info)
    genomic_coordinates_per_transcript = {}

    for exon in exon_info:
        if exon[0] not in genomic_coordinates_per_transcript:
            genomic_coordinates_per_transcript[exon[0]] = [(exon[1] - minimum, exon[2] - minimum)]
        else:
            genomic_coordinates_per_transcript[exon[0]].append((exon[1] - minimum, exon[2] - minimum))

    return genomic_coordinates_per_transcript


def get_reads_per_transcript_location(transcript_id, sqlite_path_reads):

    infile = SqliteDict(sqlite_path_reads)
    if transcript_id not in infile:
        print("No unambiguous reads support this gene")
        return None
    return infile[transcript_id]["unambig"]


def get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, filter):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.

    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_exon_coordinates = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported_transcripts, filter)
    genomic_read_dictionary = {}

    for read_file in sqlite_path_reads:
        infile = SqliteDict(read_file)

        for transcript in gene_info:
            if filter:
                if transcript[0] not in supported_transcripts:
                    continue

            if transcript[0] not in infile:
                print("No unambiguous reads support this gene")
                return None
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = genomic_exon_coordinates[transcript[0]]

            exon_junctions = string_to_integer_list(transcript[1].split(","))
            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(sequence, exons, exon_junctions)

            for length in transcript_read_dictionary:
                for location in transcript_read_dictionary[length]:
                    position = location + infile["offsets"]["fiveprime"]["offsets"][length]

                    range_counter = 0
                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_asite = position - exon[0]
                            genomic_asite = genomic_ranges[range_counter][0] + difference_between_read_position_and_exon_asite
                            if genomic_asite not in genomic_read_dictionary:
                                genomic_read_dictionary[genomic_asite] = transcript_read_dictionary[length][location]
                            else:
                                genomic_read_dictionary[genomic_asite] += transcript_read_dictionary[length][location]
                        range_counter += 1

    return genomic_read_dictionary



def get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, filter):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.

    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_exon_coordinates = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported_transcripts, filter)
    genomic_read_dictionary = {}

    for read_file in sqlite_path_reads:
        infile = SqliteDict(read_file)

        for transcript in gene_info:
            if filter:
                if transcript[0] not in supported_transcripts:
                    continue
            if transcript[0] not in infile:
                print("No unambiguous reads support this gene")
                return None
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = genomic_exon_coordinates[transcript[0]]

            exon_junctions = string_to_integer_list(transcript[1].split(","))
            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(sequence, exons, exon_junctions)

            for length in transcript_read_dictionary:
                for position in transcript_read_dictionary[length]:

                    range_counter = 0
                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_start = position - exon[0]
                            genomic_start_pos = genomic_ranges[range_counter][0] + difference_between_read_position_and_exon_start
                            if genomic_start_pos not in genomic_read_dictionary:
                                genomic_read_dictionary[genomic_start_pos] = transcript_read_dictionary[length][position]
                            else:
                                genomic_read_dictionary[genomic_start_pos] += transcript_read_dictionary[length][position]
                        range_counter += 1
    return genomic_read_dictionary


def get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, filter):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.

    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_exon_coordinates = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported_transcripts, filter)
    genomic_read_dictionary = {}

    for read_file in sqlite_path_reads:
        infile = SqliteDict(read_file)

        for transcript in gene_info:
            if filter:
                if transcript[0] not in supported_transcripts:
                    continue
            if transcript[0] not in infile:
                print("No unambiguous reads support this gene")
                return None
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = genomic_exon_coordinates[transcript[0]]

            exon_junctions = string_to_integer_list(transcript[1].split(","))
            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(sequence, exons, exon_junctions)

            for length in transcript_read_dictionary:
                for position in transcript_read_dictionary[length]:
                    range_counter = 0
                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_start = position - exon[0]
                            genomic_start_pos = genomic_ranges[range_counter][0] + difference_between_read_position_and_exon_start
                            genomic_read_range = (genomic_start_pos, genomic_start_pos + length)

                            if genomic_read_range not in genomic_read_dictionary:
                                genomic_read_dictionary[genomic_read_range] = transcript_read_dictionary[length][position]
                            else:
                                genomic_read_dictionary[genomic_read_range] += transcript_read_dictionary[length][position]
                        range_counter += 1

    return genomic_read_dictionary


def read_ranges_per_transcript(reads):
    #   Return a list of tuples of ranges which each read spans
    range_dict = {}
    for length in reads:
        for position in reads[length]:
            range_dict[(position, position + length)] = reads[length][position]


def get_exonjunction_pileup_for_transcript(transcript_id, sqlite_path_organism, sqlite_path_reads):
    # count the number of reads in the read file that span each exon-exon junction. for a given transcript
    # returns a dictionary with junctions as keys and counts as values d
    transcript_info = get_transcript_info(transcript_id, sqlite_path_organism)
    exon_junctions = string_to_integer_list(transcript_info[0][0].split(","))
    counts = {}

    for read_file in sqlite_path_reads:
        reads = get_reads_per_transcript_location(transcript_id, read_file)

        for junction in exon_junctions:
            for read_length in reads:
                for position in reads[read_length]:
                    if junction in range(position, position + read_length):
                        if junction in counts:
                            counts[junction] += reads[read_length][position]
                        else:
                            counts[junction] = reads[read_length][position]
    return counts


def get_exon_pileup_for_transcript(transcript_id, sqlite_path_organism, sqlite_path_reads):
    # count the number of reads whos p site lies within the exon sequence for a given transcript
    # returns a dictionary with sequences as keys and counts as values
    transcript_info = get_transcript_info(transcript_id, sqlite_path_organism)
    exon_junctions = string_to_integer_list(transcript_info[0][0].split(","))
    sequence = transcript_info[0][1]
    exons = exons_of_transcript(transcript_id, sqlite_path_organism)
    ranges = get_exon_coordinate_ranges(sequence, exons, exon_junctions)
    counts = {}

    for read_file in sqlite_path_reads:
        reads = get_reads_per_transcript_location(transcript_id, read_file)
        if reads is None:
            return counts

        for exon in range(len(ranges)):
            for read_length in reads:
                for position in reads[read_length]:
                    if position in range(ranges[exon][0], ranges[exon][1]+1):
                        if exons[exon] in counts:
                            counts[exons[exon]] += reads[read_length][position]
                        else:
                            counts[exons[exon]] = reads[read_length][position]
    for exon in exons:
        if exon not in counts:
            counts[exon] = 0

    return counts


def filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads):
    gene_info = get_gene_info(gene, sqlite_path_organism)
    transcripts = [i[0] for i in gene_info]
    supported = []

    for transcript in transcripts:
        pileup = get_exon_pileup_for_transcript(transcript, sqlite_path_organism, sqlite_path_reads)
        support = True

        for exon in pileup:
            if pileup[exon] == 0:
                support = False

        if support:
            supported.append(transcript)

    return supported


# def get_exon_sequences(gene, sqlite_path_organism, supported_transcripts, filter=True):
#     # returns the exon sequences for all annotated exons of a gene.
#     gene_info = get_gene_info(gene, sqlite_path_organism)
#     exon_dict = {}
#     for transcript in gene_info:
#         trans_id = transcript[0]
#         if filter and (trans_id not in supported_transcripts):
#             continue
#         exon_junctions = transcript[1].split(",")
#         sequence = transcript[2]
#         exon_junct_int = string_to_integer_list(exon_junctions)
#
#         while len(exon_junct_int) != 0:
#             exon, sequence = get_3prime_exon(exon_junct_int, sequence)
#             if exon in exon_dict:
#                 exon_dict[exon] += trans_id
#             else:
#                 exon_dict[exon] = trans_id
#             exon_junct_int.pop(-1)
#
#         if len(exon_junct_int) == 0:
#             if sequence in exon_dict:
#                 exon_dict[sequence] += trans_id
#             else:
#                 exon_dict[sequence] = trans_id
#
#         # for exon_sequence in exon_list:
#         #     if exon_sequence in exon_dict:
#         #         exon_dict[exon_sequence] += trans_id
#         #     else:
#         #         exon_dict[exon_sequence] = trans_id
#     return exon_dict

def get_exon_sequences(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # returns the exon sequences for all annotated exons of a gene.
    gene_info = get_gene_info(gene, sqlite_path_organism)
    exon_dict = {}
    for transcript in gene_info:
        trans_id = transcript[0]
        if filter and (trans_id not in supported_transcripts):
            continue
        exon_junctions = transcript[1].split(",")
        sequence = transcript[2]
        exon_junct_int = string_to_integer_list(exon_junctions)
        exon_list = []

        while len(exon_junct_int) != 0:
            exon, sequence = get_3prime_exon(exon_junct_int, sequence)
            exon_list.append(exon)
            exon_junct_int.pop(-1)

        if len(exon_junct_int) == 0:
            exon_list.append(sequence)

        # I reverse the list to fix the numbers in the genomic positioning plot to make it more intuitive.
        # exon_list = exon_list[::-1]
        for exon_sequence in exon_list:
            if exon_sequence in exon_dict:
                exon_dict[exon_sequence] += trans_id
            else:
                exon_dict[exon_sequence] = trans_id

    return exon_dict


def assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # For ease of reference each exon in the gene is given a unique integer ID
    #  retruns a dictionary with numbers as keys and sequence as values.
    exons = get_exon_sequences(gene, sqlite_path_organism, supported_transcripts, filter)
    counter = 1
    exon_dictionary = {}
    used_ids = []
    for i in exons:
        if i not in used_ids:
            exon_dictionary[counter] = i
            counter += 1
            used_ids.append(i)
    return exon_dictionary


def get_scores_per_exonjunction_for_gene(gene_name, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=True):
    # count the reads in the reads file whos p sites lite within the exon sequence
    # returns a dictionary with all unique exons in the gene as keys and counts as values
    gene_info = get_gene_info(gene_name, sqlite_path_organism)

    if filter:
        transcripts = supported_transcripts.copy()
    else:
        transcripts = [i[0] for i in gene_info]

    pileup = {}
    for transcript in transcripts:
        if transcript in pileup:
            pileup[transcript] += get_exonjunction_pileup_for_transcript(transcript, sqlite_path_organism,
                                                                         sqlite_path_reads)
        else:
            pileup[transcript] = get_exonjunction_pileup_for_transcript(transcript, sqlite_path_organism,
                                                                        sqlite_path_reads)
    return pileup


def get_scores_per_exon_for_gene(gene_name, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=True):
    # count the reads in the reads file whos p sites lite within the exon sequence
    # returns a dictionary with all unique exons in the gene as keys and counts as values
    gene_info = get_gene_info(gene_name, sqlite_path_organism)

    if filter:
        transcripts = supported_transcripts.copy()
    else:
        transcripts = [i[0] for i in gene_info]

    pileup = {}
    for transcript in transcripts:
        if transcript in pileup:
            pileup[transcript] += get_exon_pileup_for_transcript(transcript, sqlite_path_organism, sqlite_path_reads)
        else:
            pileup[transcript] = get_exon_pileup_for_transcript(transcript, sqlite_path_organism, sqlite_path_reads)

    all_exons = {}

    for transcript in pileup:
        for exon in pileup[transcript]:
            if exon in all_exons:
                all_exons[exon] += pileup[transcript][exon]
            else:
                all_exons[exon] = pileup[transcript][exon]

    return all_exons


# supported = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)
# prime = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism, supported, filter=True)
# asite = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism, supported, filter=True)
# print("asite", asite)
# print("prime", prime)
#

def get_scores_for_gene(gene_name, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=True):
    # count the reads in the reads file whos p sites lite within the exon sequence
    # returns a dictionary with all unique exons in the gene as keys and counts as values
    gene_info = get_gene_info(gene_name, sqlite_path_organism)

    if filter:
        transcripts = supported_transcripts.copy()
    else:
        transcripts = [i[0] for i in gene_info]

    exon_pileup = {}
    junct_pileup = {}
    for transcript in transcripts:
        if transcript in exon_pileup:
            exon_pileup[transcript] += get_exon_pileup_for_transcript(transcript, sqlite_path_organism, sqlite_path_reads)
        else:
            exon_pileup[transcript] = get_exon_pileup_for_transcript(transcript, sqlite_path_organism, sqlite_path_reads)

        if transcript in junct_pileup:
            junct_pileup[transcript] += get_exonjunction_pileup_for_transcript(transcript, sqlite_path_organism, sqlite_path_reads)
        else:
            junct_pileup[transcript] = get_exonjunction_pileup_for_transcript(transcript, sqlite_path_organism,
                                                                               sqlite_path_reads)
    all_exons = {}

    for transcript in exon_pileup:
        for exon in exon_pileup[transcript]:
            if exon in all_exons:
                all_exons[exon] += exon_pileup[transcript][exon]
            else:
                all_exons[exon] = exon_pileup[transcript][exon]

    exon_dict = assign_exon_numbers(gene_name, sqlite_path_organism, supported_transcripts, filter)
    scores_number = {}
    for i in all_exons:
        score_number = get_keys_by_value(exon_dict, i)[0]
        if score_number in scores_number:
            scores_number[score_number] += all_exons[i]
        else:
            scores_number[score_number] = all_exons[i]

    return scores_number, junct_pileup


def get_edges(gene_name, sqlite_path_organism, supported_transcripts, filter=True):
    # Exons that lead into eachother in a transcript are represented as directed edges in the graphs.
    # Each edge is represented as a tuple (out node, in_node). Function returns a list of edges
    gene_info = get_gene_info(gene_name, sqlite_path_organism)
    if filter:
        transcripts = supported_transcripts.copy()
    else:
        transcripts = [i[0] for i in gene_info]
    edges = []
    for transcript in transcripts:
        exon_lst = exons_of_transcript(transcript, sqlite_path_organism)
        for i in range(len(exon_lst)-1):
            edges.append((exon_lst[i], exon_lst[i+1]))
    return edges


def get_edges_scores(gene_name, sqlite_path_organism, supported_transcripts, junction_pileup, filter=True):
    # Exons that lead into eachother in a transcript are represented as directed edges in the graphs.
    # Each edge is represented as a tuple (out node, in_node). Function returns a list of edges
    gene_info = get_gene_info(gene_name, sqlite_path_organism)
    if filter:
        transcripts = supported_transcripts.copy()
    else:
        transcripts = [i[0] for i in gene_info]
    edge_scores = {}
    for transcript in transcripts:
        try:
            scores = junction_pileup[transcript]
            junct = sorted(scores.keys())
        except KeyError:
            continue
        edges = []

        exon_lst = exons_of_transcript(transcript, sqlite_path_organism)
        for i in range(len(exon_lst)-1):
            edges.append((exon_lst[i], exon_lst[i+1]))

        edge_junction = zip(edges, junct)
        for junction in edge_junction:
            edge_scores[junction[0]] = scores[junction[1]]

    return edge_scores



def get_keys_by_value(dictionary, value):
    # function to reverse look up a given dictionary returns a list of the keys with that value
    keys = []
    items = dictionary.items()
    for item in items:
        if item[1] == value:
            keys.append(item[0])
    return keys


def scores_per_exon_number(gene, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=True):
    # converts the scores dictionary from sequences as keys to numbers as keys
    scores = get_scores_per_exon_for_gene(gene, sqlite_path_organism, sqlite_path_reads, supported_transcripts, filter=filter)
    exon_dict = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter)
    scores_number = {}
    for i in scores:
        score_number = get_keys_by_value(exon_dict, i)[0]
        if score_number in scores_number:
            scores_number[score_number] += scores[i]
        else:
            scores_number[score_number] = scores[i]

    return scores_number


def get_edges_numbers(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # converts the edges to use the number ids, Returns a list of tuples (out node, in node)
    edges = get_edges(gene, sqlite_path_organism, supported_transcripts, filter)
    exon_dictionary = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter)
    number_edges = []

    for i in edges:
        number_edge = (get_keys_by_value(exon_dictionary, i[0])[0], get_keys_by_value(exon_dictionary, i[1])[0])
        number_edges.append(number_edge)
    return number_edges


def get_edges_scores_numbers(gene, sqlite_path_organism, supported_transcripts, junction_pileup, filter=True):
    # converts the edges scores to use the number ids, Returns a dictionary with tuples (out node, in node) as keys
    edges_scores = get_edges_scores(gene, sqlite_path_organism, supported_transcripts, junction_pileup, filter=True)
    exon_dictionary = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter)
    number_edge_scores = {}
    for i in edges_scores:
        number_edge = (get_keys_by_value(exon_dictionary, i[0])[0], get_keys_by_value(exon_dictionary, i[1])[0])
        number_edge_scores[number_edge] = edges_scores[i]

    return number_edge_scores


def get_path_starts(graph):
    # Get the exons that start each path through the graph. 5' exons
    edges = graph.edges()
    start_exons = list(graph.nodes())
    for edge in edges:
        if edge[1] in start_exons:
            start_exons.pop(start_exons.index(edge[1]))

    return start_exons


def principal_path_in_numbers(gene, sqlite_path_organism, supported_transcripts, filter):
    # take the principal path in exons and return the same using number ids for inclusion in the graph
    exons = exons_of_transcript(get_genes_principal(gene, sqlite_path_organism), sqlite_path_organism)
    exon_dictionary = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter)

    number_path_of_principal = []

    for item in exons:
        number_path_of_principal.append(get_keys_by_value(exon_dictionary, item)[0])

    return number_path_of_principal


def get_exon_positions(exon_info, supported_transcripts, filter=True):
    # returns the exon sequence as keys and the genomic position of the first nucleotide of the sequence as values
    transcript_positions = {}
    transcripts = []
    for exon in exon_info:
        if exon[0] not in transcripts:
            transcripts.append(exon[0])

    for i in range(len(exon_info)):
        if exon_info[i][0] not in transcript_positions:
            transcript_positions[exon_info[i][0]] = [exon_info[i][1]]
        else:
            transcript_positions[exon_info[i][0]].append(exon_info[i][1])

    exon_positions = {}

    for transcript in transcript_positions:
        exons = exons_of_transcript(transcript, sqlite_path_organism)
        positions = transcript_positions[transcript]
        for i in range(len(exons)):
            if exons[i] not in exon_positions:
                exon_positions[exons[i]] = positions[i]
            else:
                exon_positions[exons[i]] = max(positions[i], exon_positions[exons[i]])
    return exon_positions


def exon_positions_with_numbers(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # returns a dictionary of exon numbers as designated with "assign_exon_numbers()"
    # values are start and stop positions (position start, start + length)
    exon_info = get_exon_info(gene, sqlite_path_organism, supported_transcripts, filter)
    exon_positions = get_exon_positions(exon_info, supported_transcripts, filter)
    exon_dict = assign_exon_numbers(gene, sqlite_path_organism, supported_transcripts, filter)
    min_start_position = get_max_min_exon_genomic_positions(exon_info)[0]
    positions_number = {}

    for i in exon_positions:
        position_number = get_keys_by_value(exon_dict, i)[0]
        position = exon_positions[i] - min_start_position

        if position_number not in positions_number:
            positions_number[position_number] = (position, position + len(i) - 1)
        else:
            positions_number[position_number] = (max(position, positions_number[position_number][0]), 
                                                 max(position, positions_number[position_number][0]) 
                                                 + max(len(i) - 1, positions_number[position_number][1]))

    return positions_number


def get_y_values(positions_number):
    # to produce a clear layout of the networkx plot ensure that nodes dont overlap.
    # if overlap would occur based on the genomic coordinates (x axis) elevate the next
    # node one point on the y axis

    start_positions = []
    end_positions = []
    for i in positions_number:
        start_positions.append(positions_number[i][0])
        end_positions.append(positions_number[i][1])

    min_position = min(start_positions)
    max_position = max(end_positions)
    number_x_y = {}
    for i in range(min_position, max_position):
        exons_in_range = 0
        for exon in positions_number:
            if i in range(positions_number[exon][0]-100, positions_number[exon][1]+100):
                exons_in_range += 1
                if exon not in number_x_y:
                    number_x_y[exon] = [positions_number[exon][0], exons_in_range]

                else:
                    if number_x_y[exon][1] < exons_in_range:
                        number_x_y[exon][1] = exons_in_range

    return number_x_y


def relabel_nodes_L_R(number_x_y):
    relabel = {}
    node = 1
    for i in number_x_y:
        relabel[i] = node
        node += 1
    return relabel


def get_transcript_modules(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # modules are either unique or shared regions across the set of transcripts.
    # the number of exons transcribed from each individual genomic loci is recorded
    # a module is defined as a change in the number of exons transcribed from that region 
    positions_number = exon_positions_with_numbers(gene, sqlite_path_organism, supported_transcripts, filter)
    start_positions = []
    end_positions = []
    for i in positions_number:
        start_positions.append(positions_number[i][0])
        end_positions.append(positions_number[i][1])

    min_position = min(start_positions)
    max_position = max(end_positions)
    number_supported = []
    modules = {}
    module_number = 0
    module_start = min_position
    prev_count = 0

    for i in range(min_position, max_position):
        exons_in_range = 0
        exons_in_range_list = []
        for exon in positions_number:
            if i in range(positions_number[exon][0], positions_number[exon][1]):
                exons_in_range += 1
                exons_in_range_list.append(exon)

        if prev_count != exons_in_range:
            if module_number not in modules:
                modules[module_number] = (module_start, i, exons_in_range_list)
            module_number += 1
            module_start = i
        prev_count = exons_in_range

        number_supported.append(exons_in_range)

    if modules[0] == (0, 0, 0):
        modules.pop(0)

    return modules


def get_principal_modules(gene, sqlite_path_organism, modules, supported_transcripts, filter=True):
    principal_path = principal_path_in_numbers(gene, sqlite_path_organism, supported_transcripts, filter)
    modules_in_principal = []
    for module in modules:
        for exon in modules[module][2]:
            if exon in principal_path:
                modules_in_principal.append(module)

    return modules_in_principal


def get_module_pileup(modules, genomic_read_dictionary):
    # count the reads that support each module. if no reads are counted then assign a score of zero.
    module_pileup = {}
    for position in genomic_read_dictionary:
        for module in modules:
            if position in range(modules[module][0], modules[module][1]):

                if module not in module_pileup:
                    module_pileup[module] = genomic_read_dictionary[position]
                else:
                    module_pileup[module] += genomic_read_dictionary[position]

    for i, value in enumerate(modules):
        if i not in module_pileup:
            module_pileup[i] = 0
    return module_pileup


def supertranscript_positions(modules):
    # convert the modules dictionary into a list of positions of each module.
    positions = []
    for i in modules:
        positions.append((modules[i][0], modules[i][1] - modules[i][0]))

    return positions


def supertranscript_widths_starts(modules):
    # convert the modules dictionary into a list of start positions and lengths of each module.
    widths = []
    starts = []
    for i in modules:
        widths.append(modules[i][1] - modules[i][0])
        starts.append(modules[i][0])
    return widths, starts


def supertranscript_colors(scores):
    # convert the scores to a decimal between 0-1. The max being 1 and 0 values being 0.
    colors = []
    for i in range(len(scores)):
        colors.append(scores[i])
    maximum = max(colors)

    transformed_colors = []
    for i in colors:
        transformed_colors.append(i/maximum)
    return transformed_colors


def supertranscript_edgecolors(modules, principal_modules):
    # Function to provide edge color list where principal isoform path is coloured in red.
    edgecolors = []
    for i in modules:
        if i in principal_modules:
            edgecolors.append("red")
        else:
            edgecolors.append("black")

    return edgecolors


def plot_supertranscript(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # plot the super transcript using a broken bar chart with each module touching. scores/colors not supported

    modules = get_transcript_modules(gene, sqlite_path_organism, supported_transcripts, filter)
    principal_modules = get_principal_modules(gene, sqlite_path_organism, supported_transcripts, filter)
    fig, ax = plt.subplots()
    positions = supertranscript_positions(modules)
    edgecolors = supertranscript_edgecolors(modules, principal_modules)

    ax.broken_barh(positions, (0, 1), facecolors='blue', edgecolors=edgecolors)
    ax.set_ylim(0, 15)
    ax.set_xlim(0, positions[-1][0] + positions[-1][1])
    plt.show()


def plot_supertranscript_barh(gene, sqlite_path_organism, sqlite_path_reads, filter=True):
    # horizontal bar chart method of plotting a superTranscript. A modules support is
    # reflected in the color concentration.
    supported_transcripts = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)
    modules = get_transcript_modules(gene, sqlite_path_organism, supported_transcripts)
    genomic_read_dictionary = get_reads_per_genomic_location(gene, sqlite_path_reads, sqlite_path_organism,
                                                             supported_transcripts, filter)
    scores = get_module_pileup(modules, genomic_read_dictionary)
    colors = supertranscript_colors(scores)
    widths, starts = supertranscript_widths_starts(modules)
    fig, ax = plt.subplots(figsize=(10, 1.5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(False)
    c = plt.get_cmap('Blues')(colors)

    principal_modules = get_principal_modules(gene, sqlite_path_organism)
    edgecolors = supertranscript_edgecolors(modules, principal_modules)

    ax.barh(1, widths, left=starts, height=0.5, edgecolor=edgecolors, color=c)
    
    plt.show()

def get_all_paths(graph):
    # function to princt out all valid paths through a graph. Which should be
    # all valid transcripts as per the provided annotations
    for i in range(1, len(graph.nodes)+1):
        print('source: ', i)
        for j in range(1, len(graph.nodes)+1):
            print('traget: ', j)
            print(list(nx.all_simple_paths(graph, source=i, target=j)))


def nx_plot_spring(gene, sqlite_path_organism, sqlite_path_reads, filter=filter):
    graph = nx.DiGraph()
    edge_numbers = get_edges_numbers(gene, sqlite_path_organism, sqlite_path_reads, filter)
    graph.add_edges_from(edge_numbers, length=1)
    handles = []
    exon_dict = assign_exon_numbers(gene, sqlite_path_organism, sqlite_path_reads, filter)
    number_of_exons = len(exon_dict.keys())

    scores = scores_per_exon_number(gene, sqlite_path_organism, sqlite_path_reads, filter=filter)
    for i in range(number_of_exons+1):
        if i not in scores:
            scores[i] = 0
    for i in exon_dict:
        handles.append(mpatches.Patch(label=str(i) + " : "
                                      + ": " + exon_dict[i]))
        # nx.draw_networkx_labels(graph, pos) # Including nodes allowing font changes
        # Causes errors with framing the plot.
    connected = list(nx.weakly_connected_components(graph))

    color_lookup = scores
    low, *_, high = sorted(color_lookup.values())
    norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)

    gs = mpl.gridspec.GridSpec(len(connected) + 3, 1)
    plt.gca().set_position(gs[1].get_position(plt.gcf()))
    plt.gca().set_subplotspec(gs[1])

    # rcParams['figure.figsize'] = 12, 7 #full screen

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

    # G = nx.Graph()
    # plt.subplot(gs[len(connected)])
    # nx.draw(G)
    # plt.subplot(gs[len(connected)+1])
    # nx.draw(G)
    # plt.subplot(gs[len(connected)+2])
    # nx.draw(G)
    # plt.legend(loc = "lower left", prop= {"size": 4}, handles = handles)
    title = 'Kamada Kawai Layout Network of ' + gene + ' Exons'
    plt.figtext(0.0, 0.95, title)
    plt.margins(x=0, y=2)
    plt.show()
