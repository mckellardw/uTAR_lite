
import pdb
########################################################################################################
# Configfile
########################################################################################################
configfile:'STARsolo_config.yaml'
########################################################################################################
# Directories and locations
########################################################################################################
DATADIR = config['DATADIR']
TMPDIR = config['TMPDIR']

########################################################################################################
# Variables and references
########################################################################################################
GENOME_FASTA = config['GENOME_FASTA']
GENES_GTF = config['GENES_GTF']
MERGEBP = str(config['MERGEBP'])
THRESH = str(config['THRESH'])
CORES = config['CORES']

########################################################################################################
# Executables
########################################################################################################
PICARD = config['PICARD']
DROPSEQ_EXEC = config['DROPSEQ']
GTFTOGENEPRED = config['GTFTOGENEPRED']
SAMTOOLS_EXEC = config["SAMTOOLS_EXEC"]
UMITOOLS_EXEC = config["UMITOOLS_EXEC"]
FEATURECOUNTS_EXEC = config["FEATURECOUNTS_EXEC"]

########################################################################################################
rule all:
	input: 
		expand(
			'{DATADIR}/{sample}/TAR/uTAR.{FORMAT}', 
			DATADIR=config['DATADIR'], 
			sample=config['Samples'],
			FORMAT=["mtx.gz"] # ,"h5"
		)


#####################################################################################
# Set up before running HMM
#####################################################################################

# Make output directory for TAR results in STARsolo directory
rule makeOutDir:
	input:
		'{DATADIR}/{sample}',
	output:
		directory('{DATADIR}/{sample}/TAR')
	shell:
		"""
		mkdir {output}
		"""

include: "rules/1_hmm.smk"
include: "rules/2_mat.smk"
