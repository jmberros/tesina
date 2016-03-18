import pandas as pd


GENOME_STATS_FILENAME = "~/tesina/files/GRCh38.p6_genome_stats.csv"
GENOME_REGIONS_FILENAME = "~/tesina/files/GRCh39.p6_genomic_regions_definitions.txt"

def create_genome_df():
    chrom_lengths = pd.read_csv(GENOME_STATS_FILENAME, index_col="chr",
                                usecols=["chr", "total length"])
    chrom_lengths.drop(["X", "Y"], inplace=True)

    regions = pd.read_csv(GENOME_REGIONS_FILENAME, sep="\t",
                          index_col="#region_name",
                          usecols=["#region_name", "start", "stop"])
    centromers = regions.loc[[i for i in regions.index if "CEN" in i]]
    centromers["chromosome"] = [i.replace("CEN", "") for i in centromers.index]
    centromers = centromers.set_index("chromosome").drop(["X", "Y"])
    centromers.index.name = "chr"

    genome = centromers.join(chrom_lengths)
    genome.columns = ["centromere_start", "centromere_end", "chr_length"]
    genome.index = [int(i) for i in genome.index]

    return genome

