import pandas as pd
from os.path import join, expanduser


class Genome:
    BASE_DIR = expanduser("~/tesina/files")
    GENOME_STATS_FILENAME = "GRCh38.p6_genome_stats.csv"
    GENOME_REGIONS_FILENAME = "GRCh39.p6_genomic_regions_definitions.txt"

    @classmethod
    def regions(cls):
        filename = join(cls.BASE_DIR, cls.GENOME_STATS_FILENAME)
        chrom_lengths = pd.read_csv(filename, index_col="chr",
                                    usecols=["chr", "total length"])
        chrom_lengths.drop(["X", "Y"], inplace=True)

        filename = join(cls.BASE_DIR, cls.GENOME_REGIONS_FILENAME)
        regions = pd.read_csv(filename, sep="\t", index_col="#region_name",
                              usecols=["#region_name", "start", "stop"])
        centromers = regions.loc[[i for i in regions.index if "CEN" in i]]
        centromers["chromosome"] = [i.replace("CEN", "")
                                    for i in centromers.index]
        centromers = centromers.set_index("chromosome").drop(["X", "Y"])
        centromers.index.name = "chr"

        genome = centromers.join(chrom_lengths)
        genome.columns = ["centromere_start", "centromere_end", "chr_length"]
        genome.index = genome.index.astype(int)

        return genome

