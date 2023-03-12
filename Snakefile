#########################################################################
#    uTAR Snakemake pipeline - **starting from STARsolo output**
#    Copyright (C) 2020 Michael Wang
#		(modified by David McKellar)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##########################################################################

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

########################################################################################################
rule all:
############## original default call here
	input: expand('{DATADIR}/{sample}/TAR/TAR_diff_uTAR_Features_Labeled.txt', DATADIR=config['DATADIR'], sample=config['Samples'])

rule getMats:
	input: expand('{DATADIR}/{sample}/TAR/getMats.txt', DATADIR=config['DATADIR'], sample=config['Samples'])

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

########################################################################################################
# Post-alignment
########################################################################################################

include: "rules/1_hmm.smk"
include: "rules/2_mat.smk"



# from diff uTARs, extract fasta region to blast
rule getDiffRegionsToBlast:
	input:
		diffFeatures='{DATADIR}/{sample}/TAR/diff_Expressed_Features.txt',
		bamFile='{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bam',
		bamIndex='{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bam.bai'
	output:
		fastaRegions='{DATADIR}/{sample}/TAR/TAR_seqsToBlast.txt',
		uTARDiffFeatures='{DATADIR}/{sample}/TAR/TAR_diff_uTAR_Features.txt'
	shell:
		"""
		Rscript scripts/getFastasForBlast.R {input.diffFeatures} {input.bamFile}
		"""

# label uTARs based on best blast results
rule labelDiffuTARs:
	input:
		DIFF_FEAT='{DATADIR}/{sample}/TAR/TAR_diff_uTAR_Features.txt',
		BLAST_RESULT='{DATADIR}/{sample}/TAR/TAR_blastResults.txt'
	output:
		uTARSummary='{DATADIR}/{sample}/TAR/TAR_diff_uTAR_Features_Labeled.txt'
	shell:
		"""
		Rscript scripts/examineBlastResults.R {input.DIFF_FEAT} {input.BLAST_RESULT}
		"""

rule getMatsSteps:
	input:
		COUNTS='{DATADIR}/{sample}/TAR/TAR_expression_matrix_withDir.txt.gz',
		MTX_DIR=directory('{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix'),
		MTX_BC='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/barcodes.tsv.gz',
		MTX_FEAT='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/features.tsv.gz',
		MTX_MTX='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/matrix.mtx.gz',
		# hmm2='{path}/{sample}_TAR_expression_matrix_noDir.txt.gz',
		bamFile='{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bam',
		bamIndex='{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bai'
	output: 
		'{DATADIR}/{sample}/TAR/getMats.txt'
	shell:
		"""echo "Expression matrices ready" > {output}"""
