import _sqlite3
import sqlitedict
import datetime

from tripsSplice import get_reads_per_genomic_location_fiveprime
from tripsSplice import get_reads_per_genomic_location_asite
from tripsSplice import filter_unsupported_transcripts
from tripsSplice import genomic_exon_coordinate_ranges
from tripsSplice import get_read_ranges_genomic_location
from tripsSplice import get_gene_info


def exons_unique_regions(genomic_exon_coordinates):
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

                if (unique_exon[0] in range(coordinates[0], coordinates[1])) and (
                        unique_exon[1] in range(coordinates[0], coordinates[1])):
                    continue

                if (unique_exon[0] in range(coordinates[0], coordinates[1] + 1)):
                    unique_exon = (coordinates[1], unique_exon[1])

                elif (unique_exon[1] in range(coordinates[0], coordinates[1] + 1)):
                    unique_exon = (unique_exon[0], coordinates[0])

                if (unique_exon == coordinates):
                    unique_exon = (exon[1], exon[1])

            unique_regions[transcript].append(unique_exon)

    return unique_regions


def count_readranges_supporting_exons_per_transcript(exons_unique_regions, genomic_read_ranges):
    exons_counts = {}
    for read in genomic_read_ranges:
        for transcript in exons_unique_regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [0 for i in range(len(exons_unique_regions[transcript]))]

            exon_num = 0
            for exon in exons_unique_regions[transcript]:
                if exon[0] == exon[1]:
                    continue
                if (read[0] in range(exon[0], exon[1] + 1)) or (read[1] in range(exon[0], exon[1] + 1)):
                    exons_counts[transcript][exon_num] += 1
                exon_num += 1

    return exons_counts


def count_read_supporting_exons_per_transcript(exons_unique_regions, genomic_read_positions):
    exons_counts = {}
    for read in genomic_read_positions:
        for transcript in exons_unique_regions:
            if transcript not in exons_counts:
                exons_counts[transcript] = [0 for i in range(len(exons_unique_regions[transcript]))]

            exon_num = 0
            for exon in exons_unique_regions[transcript]:
                if exon[0] == exon[1]:
                    continue
                if (read[0] in range(exon[0], exon[1] + 1)) or (read[1] in range(exon[0], exon[1] + 1)):
                    exons_counts[transcript][exon_num] += 1
                exon_num += 1

    return exons_counts


def main(gene, sqlite_path_organism, sqlite_path_reads, exclude=True, count_type="range"):
    if exclude:
        supported = filter_unsupported_transcripts(gene, sqlite_path_organism, sqlite_path_reads)
    else:
        gene_info = get_gene_info(gene, sqlite_path_organism)
        supported = [i[0] for i in gene_info]
        print(supported)

    genomic_exon_coordinates = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported)
    unique_regions = exons_unique_regions(genomic_exon_coordinates)

    try:
        ["range", "fiveprime", "asite"].index(count_type)
    except ValueError:
        raise print("The count type must be one of 'range', 'fiveprime' or 'asite'. "
                    "count_type refers to the part of the read that is used in the feature counting process")

    if count_type == "range":
        genomic_read_ranges = get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported,
                                                           filter=exclude)
        counts = count_readranges_supporting_exons_per_transcript(unique_regions, genomic_read_ranges)
    if count_type == "fiveprime":
        genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_exons_per_transcript(exons_unique_regions, genomic_read_positions)
    if count_type == "asite":
        genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_exons_per_transcript(exons_unique_regions, genomic_read_positions)

    sum_of_exon_counts = {}
    maximum_sum = 0

    for transcript in counts:
        transcript_sum = sum(counts[transcript])
        sum_of_exon_counts[transcript] = transcript_sum

        if transcript_sum > maximum_sum:
            maximum_sum = transcript_sum
            maximum_transcript = transcript

    print(f"The transcript with most uniquely mapped reads is {maximum_transcript} with a score of {maximum_sum}" )
    # print(sum_of_exon_counts)

if __name__ == "__main__":

    gene = "phpt1"
    sqlite_path_organism = "homo_sapiens.v2.sqlite"
    sqlite_path_reads = ["SRR2433794.sqlite"]

    main(gene, sqlite_path_organism, sqlite_path_reads, exclude=True, count_type="asite")



