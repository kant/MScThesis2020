from orfQuant import orfQuant_signal
from orfQuant import orfQuant_OPM
import time
from tripsCountpy2 import count_read_supporting_regions_per_transcript
from tripsSplicepy2 import genomic_exon_coordinate_ranges
from tripsSplicepy2 import get_protein_coding_transcript_ids
from tripsSplicepy2 import get_reads_per_genomic_location_asite

from tripsCountpy2 import ribo_seq_read_counting
from tripsCountpy2 import ribo_seq_read_counting_raw


if __name__ == "__main__":
    start= time.time()
    gene = "phpt1"
    sqlite_path_organism = "homo_sapiens.v2.sqlite"
    sqlite_path_reads = ["SRR2433794.sqlite"]
    coding = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                  coding, filter=True)
    # genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism, coding, filter=True)
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, coding)

    counts = count_read_supporting_regions_per_transcript(exons, genomic_read_positions)

    orfQuant_res = orfQuant_OPM(gene, sqlite_path_organism, sqlite_path_reads, coding, counts, filter=True)
    end = time.time()
    print("ORFquant OPM time: " + str(end - start))
    
    start= time.time()
    gene = "phpt1"
    sqlite_path_organism = "homo_sapiens.v2.sqlite"
    sqlite_path_reads = ["SRR2433794.sqlite"]
    coding = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    genomic_read_positions = get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism,
                                                                  coding, filter=True)
    # genomic_read_positions = get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism, coding, filter=True)
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, coding)

    counts = count_read_supporting_regions_per_transcript(exons, genomic_read_positions)

    orfQuant_res = orfQuant_signal(gene, sqlite_path_organism, sqlite_path_reads, coding, counts, filter=True)
    end = time.time()
    print("ORFquant signal time: " + str(end - start))

    start= time.time()
    gene = "phpt1"
    sqlite_path_organism = "homo_sapiens.v2.sqlite"
    sqlite_path_reads = ["SRR2433794.sqlite"]

    coverage_res = ribo_seq_read_counting(gene, sqlite_path_organism, sqlite_path_reads, count_type="range", unique=True)
    end = time.time()
    print("length normalized time: " + str(end - start))

    start = time.time()
    gene = "phpt1"
    sqlite_path_organism = "homo_sapiens.v2.sqlite"
    sqlite_path_reads = ["SRR2433794.sqlite"]

    coverage_res = ribo_seq_read_counting_raw(gene, sqlite_path_organism, sqlite_path_reads, count_type="range",
                                          unique=False)
    end = time.time()
    print("raw counts time: " + str(end - start))