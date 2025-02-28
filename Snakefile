import pdb


########################################################################################################
# Configfile
########################################################################################################
configfile: "STARsolo_config.yaml"


########################################################################################################
# Directories and locations
########################################################################################################
DATADIR = config["DATADIR"]
TMPDIR = config["TMPDIR"]

########################################################################################################
# Variables and references
########################################################################################################
GENOME_FASTA = config["GENOME_FASTA"]
GENES_GTF = config["GENES_GTF"]
MERGEBP = str(config["MERGEBP"])
THRESH = str(config["THRESH"])
CORES = config["CORES"]


########################################################################################################
rule all:
    input:
        expand(
            "{DATADIR}/{sample}/TAR/uTAR.{FORMAT}",
            DATADIR=config["DATADIR"],
            sample=config["Samples"],
            FORMAT=["mtx.gz"],
        ),


#####################################################################################
# Set up before running HMM
#####################################################################################


# Make output directory for TAR results in STARsolo directory
rule makeOutDir:
    input:
        "{DATADIR}/{sample}",
    output:
        directory("{DATADIR}/{sample}/TAR"),
    shell:
        """
        mkdir {output}
        """


# include: "rules/0_dedup.smk"
include: "rules/1_hmm.smk"
include: "rules/2_mat.smk"


#### Utility rules ############################
# Index .bam file
rule utils_index_BAM:
    input:
        BAM="{BAM}",
    output:
        BAI="{BAM}.bai",
    # wildcard_constraints:
    #     BAM=".*\.(bam)$"
    # resources:
    threads: config["CORES"]
    shell:
        """
        samtools index -@ {threads} {input.BAM}
        """
