from tripsSplicepy2 import genomic_exon_coordinate_ranges
from tripsSplicepy2 import get_protein_coding_transcript_ids
from tripsSplicepy2 import genomic_orf_coordinate_ranges
from tripsSplicepy2 import genomic_junction_positions
from tripsSplicepy2 import get_reads_per_transcript_location


from tripsSplicepy2 import get_reads_per_genomic_location_fiveprime
from tripsSplicepy2 import get_reads_per_genomic_location_asite
from tripsSplicepy2 import get_read_ranges_genomic_location
from tripsSplicepy2 import genomic_junction_scores
from tripsSplicepy2 import get_start_stop_codon_positions

from tripsCountpy2 import count_readranges_supporting_exons_per_transcript
from tripsCountpy2 import count_read_supporting_regions_per_transcript
from tripsCountpy2 import rank_based_on_dict_values
import time



def classify_regions_pos_neg(regions, supported, exclude=True, count_type="range"):
    # classify the genomic coordinate ranges into pos or neg. pos if 1 or more read supports it. Returns a nested dictionary
    # with two nests. 1st on transcripts as keys, 2nd with pos or neg as keys

    if count_type == "range":
        genomic_read_ranges = get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported,
                                                           filter=exclude)
        counts = count_readranges_supporting_exons_per_transcript(regions, genomic_read_ranges)

    if count_type == "fiveprime":
        genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(regions, genomic_read_positions)

    if count_type == "asite":
        genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                          supported, filter=exclude)
        counts = count_read_supporting_regions_per_transcript(regions, genomic_read_positions)


    classified = {}
    for transcript in counts:
        if transcript not in classified:
            classified[transcript] = {'pos':[], 'neg':[]}

        for i in zip(regions[transcript], counts[transcript]):
            if i[1] > 0:
                classified[transcript]['pos'].append(i[0])
            else:
                classified[transcript]['neg'].append(i[0])
    return classified


def classify_regions_shared_unique(regions):
    classified = {}
    for tran1 in regions:
        if tran1 not in classified:
            classified[tran1] = {"unique":[], "shared":[]}
        for region in regions[tran1]:
            counter = 0
            for tran2 in regions:
                counter += regions[tran2].count(region)

            if counter > 1:
                classified[tran1]["shared"].append(region)
            else:
                classified[tran1]["unique"].append(region)
    return classified


def classify_junctions_pos_neg(junction_scores):
    classified = {}

    for transcript in junction_scores:
        if transcript not in classified:
            classified[transcript] = {"pos":[], "neg":[]}

        for junction in junction_scores[transcript]:
            if junction_scores[transcript][junction] > 0:

                if junction not in classified[transcript]:
                    classified[transcript]['pos'].append(junction)

            else:
                if junction not in classified[transcript]:
                    classified[transcript]['neg'].append(junction)

    return classified


def get_orf_coordinate_junctions(orf_coordinates):

    junctions = {}
    for transcript in orf_coordinates:
        if transcript not in junctions:
            junctions[transcript] = []
        for i in range(len(orf_coordinates[transcript]) - 1):
            junctions[transcript].append((orf_coordinates[transcript][i][1], orf_coordinates[transcript][i + 1][0] - 1))

    return junctions


def features_per_orf(orf_coordinates, orf_junctions, exons, junctions):
    features_in_orf = {}

    for transcript in orf_coordinates:
        if transcript not in features_in_orf:
            features_in_orf[transcript] = {"junctions":[], "exons":[]}

            for region in orf_coordinates[transcript]:
                for exon in exons[transcript]:
                    if (exon == region):
                        features_in_orf[transcript]["exons"].append(exon)

                    elif (region[1] in range(exon[0], exon[1])) and region[0] == exon[0]:
                        features_in_orf[transcript]["exons"].append(exon)
                    elif (region[0] in range(exon[0], exon[1])) and region[1] == exon[1]:
                        features_in_orf[transcript]["exons"].append(exon)

            for junction in orf_junctions[transcript]:
                for exon_junction in junctions[transcript]:
                    if (junction == exon_junction):
                        features_in_orf[transcript]["junctions"].append(junction)
    return features_in_orf


def region_coverage(regions, counts):
    region_counts = {}
    for transcript in regions:
        region_counts_transcript = zip(regions[transcript], counts[transcript])
        if transcript not in region_counts:
            region_counts[transcript] = region_counts_transcript

    coverage_region_transcript = {}
    for transcript in region_counts:
        if transcript not in coverage_region_transcript:
            coverage_region_transcript[transcript] = {}

        for item in region_counts[transcript]:
            length = item[0][1] - item[0][0]
            coverage = round(float(item[1])/float(length), 4)
            coverage_region_transcript[transcript][item[0]] = coverage

    return coverage_region_transcript


def coverage_junction_transcript(junction_scores):
    #function to calculate the normalised coverage of junction features. returns a nested dict.
    # 60 = 30 + 30 which is the sum of two riboseq read lengths.
    coverage_junction = {}
    for transcript in junction_scores:
        if transcript not in coverage_junction:
            coverage_junction[transcript] = {}
        for junction in junction_scores[transcript]:
            coverage = float(junction_scores[transcript][junction]) / float(60)
            coverage_junction[transcript][junction] = coverage

    return coverage_junction


def average_unique_coverage(unique_shared_exons, unique_shared_junctions, coverage_junction, coverage_exons):
    average_features = {}
    for transcript in unique_shared_exons:
        sum = 0
        count = 0
        for item in unique_shared_exons[transcript]['unique']:
            sum += coverage_exons[transcript][item]
            count += 1
        for junction in unique_shared_junctions[transcript]['unique']:
            sum += coverage_junction[transcript][junction]
            count += 1

        if transcript not in average_features:
            average_features[transcript] = float(sum)/float(count)

    return average_features


def all_feature_average(coverage_junction, coverage_exons):
    average_features = {}

    for transcript in coverage_exons:
        sum = 0
        count = 0
        for exon in coverage_exons[transcript]:
            sum += coverage_exons[transcript][exon]
            count += 1

        for junction in coverage_junction[transcript]:
            sum += coverage_junction[transcript][junction]
            count += 1


        if transcript not in average_features:
            average_features[transcript] = float(sum) / float(count)
    return average_features


def cORF_ratio(average_unique, average_all):
    cORF ={}
    for transcript in average_all:
        cORF[transcript] = average_unique[transcript] / average_all[transcript]

    return cORF


def adjusted_coverage_for_non_unique_orfs(transcript, coverage_exons, coverage_junctions, average_coverage_all):
    adjusted_coverage_exons = {transcript: {}}
    for exon in coverage_exons[transcript]:
        signal = 0
        for trans in coverage_exons:
            if exon in coverage_exons[trans]:
                signal += coverage_exons[trans][exon]

        if exon not in adjusted_coverage_exons[transcript]:
            adjusted_coverage_exons[transcript][exon] = coverage_exons[transcript][exon] - (
                    coverage_exons[transcript][exon]  * signal)

    adjusted_coverage_junctions = {transcript: {}}
    for junction in coverage_junctions[transcript]:
        signal = 0
        for trans in coverage_junctions:
            if junction in coverage_junctions[trans]:
                signal += coverage_junctions[trans][junction]

        if junction not in adjusted_coverage_junctions[transcript]:
            adjusted_coverage_junctions[transcript][junction] = coverage_junctions[transcript][junction] - (
                        coverage_junctions[transcript][junction] * signal)
    sum = 0
    count = 0
    for region in adjusted_coverage_exons[transcript]:
        sum += adjusted_coverage_exons[transcript][region]
        count += 1

    for junction in adjusted_coverage_junctions[transcript]:
        sum += adjusted_coverage_junctions[transcript][junction]
        count += 1
    average_adjusted_coverage = float(sum)/float(count)

    cORF_transcript = average_adjusted_coverage / average_coverage_all[transcript]
    return cORF_transcript


def shared_coverage(coverage_exons, coverage_junctions):
    shared_exon_coverage = {}
    for transcript in coverage_exons:
        if transcript not in shared_exon_coverage:
            shared_exon_coverage[transcript] = {}
        for exon in coverage_exons[transcript]:
            number_ORF_over_F = 0
            for tran in coverage_exons:
                for exon2 in coverage_exons[tran]:
                    if transcript == tran:
                        continue

                    if exon == exon2:
                        number_ORF_over_F += 1
            number_ORF_over_F += 1

            shared_exon_coverage[transcript][exon] = coverage_exons[transcript][exon]/number_ORF_over_F

    shared_junction_coverage = {}
    for transcript in coverage_junctions:
        if transcript not in shared_junction_coverage:
            shared_junction_coverage[transcript] = {}
        for junction in coverage_junctions[transcript]:
            number_ORF_over_F = 0
            for tran in coverage_junctions:
                for junction2 in coverage_junctions[tran]:
                    if transcript == tran:
                        continue

                    if junction == junction2:
                        number_ORF_over_F += 1
            number_ORF_over_F += 1
            shared_junction_coverage[transcript][junction] = coverage_junctions[transcript][junction] / number_ORF_over_F

    shared_cORF = {}
    for transcript in shared_exon_coverage:
        sum = 0
        count = 0
        for exon in shared_exon_coverage[transcript]:
            sum += shared_exon_coverage[transcript][exon]
            count += 1

        for junction in shared_junction_coverage[transcript]:
            sum += shared_junction_coverage[transcript][junction]
            count += 1

        if transcript not in shared_cORF:
            shared_cORF[transcript] = float(sum)/float(count)

    return shared_cORF


def aORF(cORF, counts):
    a_sites_perORF = {}
    for transcript in cORF:
        a = []

        if transcript not in a_sites_perORF:
            a_sites_perORF[transcript] = cORF[transcript] * sum(counts[transcript])

    return a_sites_perORF


def lORF(coding, sqlite_path_organism):
    lORFs = {}

    for transcript in coding:
        start, stop = get_start_stop_codon_positions(transcript, sqlite_path_organism)
        length = (stop + 1) - start
        if transcript not in lORFs:
            lORFs[transcript] = length
    return lORFs

def ORFs_per_million(aORF, lORF):
    aORF_lORF = {}
    full_set = 0
    for transcript in aORF:
        if transcript not in aORF_lORF:
            aORF_lORF[transcript] = aORF[transcript]/lORF[transcript]
            full_set += aORF[transcript]/lORF[transcript]

    orfs_per_million = {}
    for transcript in aORF_lORF:
        if transcript not in orfs_per_million:
            orfs_per_million[transcript] = aORF_lORF[transcript] * (10*6 / full_set)
    return orfs_per_million


def pct_gene_signal_per_orf(aORF):
    total_aORF = 0
    for transcript in aORF:
        total_aORF += aORF[transcript]

    pct_ORF = {}
    for transcript in aORF:
        if transcript not in pct_ORF:
            pct_ORF[transcript] = aORF[transcript] / total_aORF

    return pct_ORF




def orfQuant_OPM(gene, sqlite_path_organism, sqlite_path_reads, supported, counts, filter=True):
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported)
    junctions = genomic_junction_positions(gene, sqlite_path_organism, supported)
    orf_coordinates = genomic_orf_coordinate_ranges(gene, sqlite_path_organism, supported)
    orf_junctions = get_orf_coordinate_junctions(orf_coordinates)
    features_per_orf(orf_coordinates, orf_junctions, exons, junctions)
    junction_scores = genomic_junction_scores(gene, sqlite_path_organism, sqlite_path_reads, supported,
                               filter=True)


    unique_shared_exons = classify_regions_shared_unique(exons)
    unique_shared_junctions = classify_regions_shared_unique(junctions)
    coverage_exons = region_coverage(exons, counts)
    coverage_junctions = coverage_junction_transcript(junction_scores)

    shared_coverage(coverage_exons, coverage_junctions)

    # pos_neg_exons = classify_regions_pos_neg(exons, supported, exclude=True, count_type="range")
    # pos_neg_junctions = classify_junctions_pos_neg(junction_scores)


    average_unique = average_unique_coverage(unique_shared_exons, unique_shared_junctions, coverage_junctions, coverage_exons)
    average_all = all_feature_average(coverage_junctions, coverage_exons)

    cORF = cORF_ratio(average_unique, average_all)

    all_shared = True
    for transcript in supported:
        if unique_shared_junctions[transcript]["unique"] == [] and unique_shared_exons[transcript]["unique"] == []:
            cORF[transcript] = adjusted_coverage_for_non_unique_orfs(transcript, coverage_exons, coverage_junctions, average_all)
        else:
            all_shared = False

    if all_shared:
        cORF = shared_coverage(coverage_exons, coverage_junctions)
    else:
        cORF = cORF_ratio(average_unique, average_all)

    adjusted_a_sites = aORF(cORF, counts)
    orf_lengths = lORF(supported, sqlite_path_organism)

    orfs_per_million = ORFs_per_million(adjusted_a_sites, orf_lengths)
    return orfs_per_million
    # gene_signal_per_orf = pct_gene_signal_per_orf(adjusted_a_sites)

    # return gene_signal_per_orf

def orfQuant_signal(gene, sqlite_path_organism, sqlite_path_reads, supported, counts, filter=True):
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported)
    junctions = genomic_junction_positions(gene, sqlite_path_organism, supported)
    orf_coordinates = genomic_orf_coordinate_ranges(gene, sqlite_path_organism, supported)
    orf_junctions = get_orf_coordinate_junctions(orf_coordinates)
    features_per_orf(orf_coordinates, orf_junctions, exons, junctions)
    junction_scores = genomic_junction_scores(gene, sqlite_path_organism, sqlite_path_reads, supported,
                               filter=True)


    unique_shared_exons = classify_regions_shared_unique(exons)
    unique_shared_junctions = classify_regions_shared_unique(junctions)
    coverage_exons = region_coverage(exons, counts)
    coverage_junctions = coverage_junction_transcript(junction_scores)

    shared_coverage(coverage_exons, coverage_junctions)

    # pos_neg_exons = classify_regions_pos_neg(exons, supported, exclude=True, count_type="range")
    # pos_neg_junctions = classify_junctions_pos_neg(junction_scores)


    average_unique = average_unique_coverage(unique_shared_exons, unique_shared_junctions, coverage_junctions, coverage_exons)
    average_all = all_feature_average(coverage_junctions, coverage_exons)

    cORF = cORF_ratio(average_unique, average_all)

    all_shared = True
    for transcript in supported:
        if unique_shared_junctions[transcript]["unique"] == [] and unique_shared_exons[transcript]["unique"] == []:
            cORF[transcript] = adjusted_coverage_for_non_unique_orfs(transcript, coverage_exons, coverage_junctions, average_all)
        else:
            all_shared = False

    if all_shared:
        cORF = shared_coverage(coverage_exons, coverage_junctions)
    else:
        cORF = cORF_ratio(average_unique, average_all)

    adjusted_a_sites = aORF(cORF, counts)
    orf_lengths = lORF(supported, sqlite_path_organism)

    gene_signal_per_orf = pct_gene_signal_per_orf(adjusted_a_sites)

    return gene_signal_per_orf

# if __name__ == "__main__":
#     start= time.time()
#     gene = "igf2"
#     sqlite_path_organism = "homo_sapiens.v2.sqlite"
#     sqlite_path_reads = ["SRR2433794.sqlite"]
#     coding = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
#     genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
#                                                                   coding, filter=True)
#     # genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism, coding, filter=True)
#     exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, coding)
#
#     counts = count_read_supporting_regions_per_transcript(exons, genomic_read_positions)
#
#     orfQuant_res = orfQuant_OPM(gene, sqlite_path_organism, sqlite_path_reads, coding, counts, filter=True)
#     end = time.time()
#     print(end - start)
#     rankings = rank_based_on_dict_values(orfQuant_res)
#     print rankings