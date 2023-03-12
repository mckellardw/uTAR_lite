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
	input: 
		expand(
			'{DATADIR}/{sample}/TAR/uTAR_matrix.mtx', 
			DATADIR=config['DATADIR'], 
			sample=config['Samples']
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
