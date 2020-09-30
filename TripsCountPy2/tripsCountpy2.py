import sqlite3

from tripsSplicepy2 import get_reads_per_genomic_location_fiveprime
from tripsSplicepy2 import get_reads_per_genomic_location_asite
from tripsSplicepy2 import get_read_ranges_genomic_location
from tripsSplicepy2 import get_scores_per_exonjunction_for_gene

from tripsSplicepy2 import filter_unsupported_transcripts
from tripsSplicepy2 import genomic_exon_coordinate_ranges
from tripsSplicepy2 import get_protein_coding_transcript_ids
from tripsSplicepy2 import genomic_orf_coordinate_ranges
from tripsSplicepy2 import get_keys_by_value

import sys
# print sys.version


def get_gene_info(gene_id, sqlite_path_organism):
    # get the transcript ID, exon junctions and transcript sequence from the specified sqlite file for a given gene ID
    gene_id = gene_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript, exon_junctions, sequence FROM transcripts WHERE gene = (:gene);', {'gene': gene_id})
    gene_info = c.fetchall()
    return gene_info


def unique_exon_junctions(genomic_exon_coordinates):
    transcript_junctions = {}
    junction_list = []
    for transcript in genomic_exon_coordinates:
        number_of_exons = len(genomic_exon_coordinates[transcript])
        for i in range(number_of_exons - 1):
            junction_list.append((genomic_exon_coordinates[transcript][i][1], genomic_exon_coordinates[transcript][i + 1][0]))

            if str(transcript) in transcript_junctions:
                transcript_junctions[str(transcript)].append((genomic_exon_coordinates[transcript][i][1], genomic_exon_coordinates[transcript][i + 1][0]))
            else:
                transcript_junctions[str(transcript)] = [(genomic_exon_coordinates[transcript][i][1], genomic_exon_coordinates[transcript][i + 1][0])]

    unique_transcript_junctions = {key: value[:] for key, value in transcript_junctions.items()}

    for junct in junction_list:
        if junction_list.count(junct) > 1:
            for transcript in unique_transcript_junctions:
                if junct in unique_transcript_junctions[transcript]:
                    ind = unique_transcript_junctions[transcript].index(junct)
                    unique_transcript_junctions[transcript].pop(ind)

    return transcript_junctions, unique_transcript_junctions


def get_unique_regions(genomic_exon_coordinates):
    # return a dictionary with ensembl ids as keys and a list of tuples of the coordinates of the regions unique to each
    # exon of that transcript. This is achieved by comparing each exon which each other exon and pairing off the shared regions
    # one exact match between exon and coordinate is okay (matching with itself) greater than this is due to the exon
    # occuring in more than one transcript

    unique_regions = {}
    exon_coordinates = []

    for transcript in genomic_exon_coordinates:
        for exon in genomic_exon_coordinates[transcript]:
            exon_coordinates.append(exon)

    for transcript in genomic_exon_coordinates:
        if transcript not in unique_regions:
            unique_regions[transcript] = []
        for exon in genomic_exon_coordinates[transcript]:
            unique_exon = exon
            exact_match_counter = 0

            for coordinates in exon_coordinates:
                if exon == coordinates:
                    exact_match_counter += 1
                    if exact_match_counter > 1:
                        unique_exon = (exon[1], exon[1])
                    continue
                if (unique_exon[0] in range(coordinates[0], coordinates[1] + 1)) and (
                        unique_exon[1] in range(coordinates[0], coordinates[1] + 1)):
                    unique_exon = (unique_exon[1], unique_exon[1])
                    continue

                if (unique_exon[0] in range(coordinates[0], coordinates[1] + 1)):
                    unique_exon = (coordinates[1], unique_exon[1])

                if (unique_exon[1] in range(coordinates[0], coordinates[1] + 1)):
                    unique_exon = (unique_exon[0], coordinates[0])

                if (unique_exon == coordinates):
                    unique_exon = (exon[1], exon[1])

            unique_regions[transcript].append(unique_exon)
    return unique_regions


def count_readranges_supporting_exons_per_transcript(regions, genomic_read_ranges):
    exons_counts = {}
    for read in genomic_read_ranges:
        for transcript in regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [0 for i in range(len(regions[transcript]))]

            exon_num = 0
            for exon in regions[transcript]:

                if exon[0] == exon[1]:
                    exon_num += 1
                    continue


                if (read[0] in range(exon[0], exon[1] + 1)) or (read[1] in range(exon[0], exon[1] + 1)):
                    # print read, exon, regions[transcript][exon_num], genomic_read_ranges[read]

                    exons_counts[transcript][exon_num] += genomic_read_ranges[read]
                exon_num += 1
            # print "regions: ", regions[transcript]
            # print "counts; ", exons_counts[transcript]

    return exons_counts


def count_read_supporting_regions_per_transcript(regions, genomic_read_positions):
    exons_counts = {}
    for read in genomic_read_positions:
        for transcript in regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [0 for i in range(len(regions[transcript]))]

            exon_num = 0
            for exon in regions[transcript]:
                if exon[0] == exon[1]:
                    exon_num += 1
                    continue
                if (read in range(exon[0], exon[1] + 1)):
                    exons_counts[transcript][exon_num] += 1
                exon_num += 1

    return exons_counts


def explain_exon_junctions(junction_scores, all_junctions, unique_junctions):
    # first this function assigns the junctions scores that are known for each junction in a transcript to those same
    # junctions but on genomic coordinate level. This is achieved by handling each transcript at a time and assigning the
    # scores based on the order they come on the genome/ transcript
    # it then returns the transcripts that are required to explain the reads that map to these unique transcript

    genomic_junction_scores = {}
    for transcript in junction_scores:
        number_of_junctions = len(junction_scores[transcript])
        score_keys = sorted(junction_scores[transcript].keys())
        if transcript not in genomic_junction_scores:
            genomic_junction_scores[transcript] = {}

        for i in range(number_of_junctions):
            if all_junctions[transcript][i] not in genomic_junction_scores[transcript]:
                genomic_junction_scores[transcript][all_junctions[transcript][i]] = junction_scores[transcript][score_keys[i]]
    required_transcripts_to_explain_reads = []
    for transcript in all_junctions:
        if len(unique_junctions[transcript]) > 0:
            for junction in unique_junctions[transcript]:
                if genomic_junction_scores[transcript][junction] > 0:
                    required_transcripts_to_explain_reads.append(transcript)

    return required_transcripts_to_explain_reads


def get_coverage_per_region(regions, counts):

    region_lengths = {}
    for transcript in regions:
        if transcript not in region_lengths:
            region_lengths[transcript] = []

        for region in regions[transcript]:
            region_lengths[transcript].append(region[1] - region[0])


    #coverage = reads/length of region
    region_coverage = {}
    for transcript in counts:
        if transcript not in region_coverage:
            region_coverage[transcript] = []

        for region, i in enumerate(counts[transcript]):
            if region_lengths[transcript][region] == 0:
                region_coverage[transcript].append(None)
            else:
                region_coverage[transcript].append(round(float(counts[transcript][region]) / float(region_lengths[transcript][region]), 4))
    return region_coverage


def average_coverage_per_transcript(region_coverage):
    average_coverage = {}
    for transcript in region_coverage:
        sum = 0
        count = 0
        if (region_coverage[transcript] == None) and all(region_coverage[transcript]):
            print "transcript {transcript} had no unique regions".format(transcript=transcript)
            if transcript not in average_coverage:
                average_coverage[transcript] = None
        else:
            for region in region_coverage[transcript]:
                if region != None:
                    count += 1
                    sum += region

        if transcript not in average_coverage:
            average_coverage[transcript] = sum/count

    return average_coverage


def rank_based_on_dict_values(transcript_dict):
    region_scores = []
    for transcript in transcript_dict:
        region_scores.append(transcript_dict[transcript])
    ordered_region_scores = sorted(region_scores, reverse=True)
    rankings = []
    for score in ordered_region_scores:
        scores_transcript = get_keys_by_value(transcript_dict, score)[0]
        rankings.append((scores_transcript, score))

    return rankings


def rna_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, exclude=True, count_type="range"):
    if exclude:
        supported = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)
    else:
        gene_info = get_gene_info(gene, sqlite_path_organism)
        supported = [i[0] for i in gene_info]

    genomic_exon_coordinates = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported)
    unique_regions = get_unique_regions(genomic_exon_coordinates)
    all_junctions, unique_junctions = unique_exon_junctions(genomic_exon_coordinates)
    junction_scores = get_scores_per_exonjunction_for_gene(gene, sqlite_path_organism, sqlite_path_reads, supported, filter=exclude)
    transcripts_to_explain_reads = explain_exon_junctions(junction_scores, all_junctions, unique_junctions)


    try:
        ["range", "fiveprime", "asite"].index(count_type)
    except ValueError:
        print("The count type must be one of 'range', 'fiveprime' or 'asite'. "
                "count_type refers to the part of the read that is used in the feature counting process")
        return "ERROR"

    if count_type == "range":
        genomic_read_ranges = get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported,
                                                           filter=exclude)
        counts = count_readranges_supporting_exons_per_transcript(unique_regions, genomic_read_ranges)

    if count_type == "fiveprime":
        genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(unique_regions, genomic_read_positions)

    if count_type == "asite":
        genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(unique_regions, genomic_read_positions)


    sum_of_exon_counts = {}
    maximum_sum = 0

    for transcript in counts:
        transcript_sum = sum(counts[transcript])
        if (transcript_sum > 0) and (transcript not in transcripts_to_explain_reads):
            transcripts_to_explain_reads.append(str(transcript))
        sum_of_exon_counts[transcript] = transcript_sum

        if transcript_sum > maximum_sum:
            maximum_sum = transcript_sum
            maximum_transcript = transcript

    print("The transcript with most uniquely mapped reads is {maximum_transcript} with a score of {maximum_sum}".format(maximum_transcript=str(maximum_transcript), maximum_sum=maximum_sum) )
    return transcripts_to_explain_reads


def ribo_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, count_type="range", unique=True):
    supported = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    exclude = True
    orf_regions = genomic_orf_coordinate_ranges(gene, sqlite_path_organism, supported)
    if unique:
        orf_regions = get_unique_regions(orf_regions)

    # unique_regions = get_unique_regions(genomic_exon_coordinates)
    # all_junctions, unique_junctions = unique_exon_junctions(genomic_exon_coordinates)
    # junction_scores = get_scores_per_exonjunction_for_gene(gene, sqlite_path_organism, sqlite_path_reads, supported, filter=exclude)
    # transcripts_to_explain_reads = explain_exon_junctions(junction_scores, all_junctions, unique_junctions)

    try:
        ["range", "fiveprime", "asite"].index(count_type)
    except ValueError:
        print("The count type must be one of 'range', 'fiveprime' or 'asite'. "
                "count_type refers to the part of the read that is used in the feature counting process")
        return "ERROR"

    if count_type == "range":
        genomic_read_ranges = get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported,
                                                           filter=exclude)
        counts = count_readranges_supporting_exons_per_transcript(orf_regions, genomic_read_ranges)

    if count_type == "fiveprime":
        genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(orf_regions, genomic_read_positions)

    if count_type == "asite":
        genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(orf_regions, genomic_read_positions)

    sum_of_exon_counts = {}
    maximum_sum = 0

    region_coverage = get_coverage_per_region(orf_regions, counts)
    coverage = average_coverage_per_transcript(region_coverage)
    return coverage
    # rankings = rank_based_on_dict_values(coverage)
    #
    # for transcript in counts:
    #     transcript_sum = sum(counts[transcript])
    #     # if (transcript_sum > 0) and (transcript not in transcripts_to_explain_reads):
    #     #     transcripts_to_explain_reads.append(str(transcript))
    #     sum_of_exon_counts[transcript] = transcript_sum
    #
    #     if transcript_sum > maximum_sum:
    #         maximum_sum = transcript_sum
    #         maximum_transcript = transcript
    # print "highest sum of counts: {maximum_transcript} with a total sum of: {maximum_sum}".format(maximum_transcript=str(maximum_transcript), maximum_sum=maximum_sum)
    # return rankings

def ribo_seq_read_counting_raw(gene, sqlite_path_organism, sqlite_path_reads, count_type="range", unique=True):
    supported = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    exclude = True
    orf_regions = genomic_orf_coordinate_ranges(gene, sqlite_path_organism, supported)
    if unique:
        orf_regions = get_unique_regions(orf_regions)

    # unique_regions = get_unique_regions(genomic_exon_coordinates)
    # all_junctions, unique_junctions = unique_exon_junctions(genomic_exon_coordinates)
    # junction_scores = get_scores_per_exonjunction_for_gene(gene, sqlite_path_organism, sqlite_path_reads, supported, filter=exclude)
    # transcripts_to_explain_reads = explain_exon_junctions(junction_scores, all_junctions, unique_junctions)

    try:
        ["range", "fiveprime", "asite"].index(count_type)
    except ValueError:
        print("The count type must be one of 'range', 'fiveprime' or 'asite'. "
                "count_type refers to the part of the read that is used in the feature counting process")
        return "ERROR"

    if count_type == "range":
        genomic_read_ranges = get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported,
                                                           filter=exclude)
        counts = count_readranges_supporting_exons_per_transcript(orf_regions, genomic_read_ranges)

    if count_type == "fiveprime":
        genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(orf_regions, genomic_read_positions)

    if count_type == "asite":
        genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(orf_regions, genomic_read_positions)

    sum_of_exon_counts = {}
    maximum_sum = 0
    print(counts)
    return counts


if __name__ == "__main__":

    gene = "phpt1"
    sqlite_path_organism = "homo_sapiens.v2.sqlite"
    sqlite_path_reads = ["SRR2433794.sqlite"]
    rankings = ribo_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, count_type="asite", unique = True)
    print rankings
    # transcripts = rna_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, exclude=False, count_type="range")
    # print "To explain all reads in the selected files must include: ", transcripts


