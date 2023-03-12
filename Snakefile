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

#FULLSCRIPTPATH=str(config['FULLSCRIPTPATH'])
########################################################################################################
# Executables
########################################################################################################
PICARD = config['PICARD']
DROPSEQ = config['DROPSEQ']
GTFTOGENEPRED = config['GTFTOGENEPRED']

########################################################################################################
rule all:
############## original default call here
	input: expand('{DATADIR}/{sample}/TAR/TAR_diff_uTAR_Features_Labeled.txt', DATADIR=config['DATADIR'], sample=config['Samples'])

rule getMats:
	input: expand('{DATADIR}/{sample}/TAR/getMats.txt', DATADIR=config['DATADIR'], sample=config['Samples'])

#####################################################################################
# Set up before running HMM
#####################################################################################

# create sequence dictionary using picard tools - only need to run once, saves into your cellranger reference
rule createGenomeDict:
	input:
		FASTA = GENOME_FASTA
	output:
		DICT = GENOME_FASTA+'.dict'
	shell:
		"""
		java -Dpicard.useLegacyParser=false -jar {PICARD} CreateSequenceDictionary \
		-R {input.FASTA} \
		-O {output.DICT}
		"""

# convert GTF to REFFlat, save in your cellranger reference
# 	not run if REFFlat file exists -  will only need to run this once for each reference genome
rule convertToRefFlat:
	input:
		GENES = GENES_GTF
	output:
		REFFLAT = GENES_GTF+".refFlat"
	shell:
		"""
		{GTFTOGENEPRED} -genePredExt -geneNameAsName2 {input.GENES} refFlat.tmp
		paste <(cut -f 12 refFlat.tmp) <(cut -f 1-10 refFlat.tmp) > {output.REFFLAT}
		rm refFlat.tmp
		"""

# Make output directory for TAR results in STARsolo directory
#@done
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

# if combined file is too big, need to downsample
#rule# downSampleBam:
#	input: "combined_bam.bam"
#	output: "combined_bam.bam"
#	shell: """samtools view -s {downSamp} -b {input} > {output}

# Run HMM
rule calcHMMbed:
	input:
		'{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam' #need absolute path
	output:
		temp('{DATADIR}/{sample}/TAR/TAR_reads.bed.gz')
	threads:
		CORES
	params:
		MEM="64G" # Change this if you are using a low-RAM machine
	shell:
		"""
		bash scripts/SingleCellHMM_MW.bash {input} {threads} {params.MEM} {MERGEBP} {THRESH} {wildcards.DATADIR}/{wildcards.sample}/TAR
		"""
# Generate annotations from HMM in refFlat format
#		Returns annotations with and without considering strandedness/direction
rule calcHMMrefFlat:
	input:
		BEDGZ='{DATADIR}/{sample}/TAR/TAR_reads.bed.gz',
		REFFLAT=GENES_GTF+".refFlat"
	output:
		# NODIR = '{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.noDir.refFlat',
		WITHDIR ='{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.withDir.refFlat'
	threads:
		CORES
	shell:
		"""
		Rscript scripts/generate_refFlat_script_both.R {input.REFFLAT} {input.BEDGZ} {threads}
		"""

# Convert stranded annotations to gtf format
rule HMM_refFlat_to_gtf_WITHDIR:
	input:
		WITHDIR ='{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.withDir.refFlat'
	output:
		WITHDIR ='{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.withDir.refFlat.gtf'
	shell:
		"""
		Rscript scripts/convertRefFlatToGTF.R {input.WITHDIR}
		"""

########################################################################################################
## Tagging and labeling
##	source: https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html
########################################################################################################

# Label .bam file with each HMM feature
#TODO: add featureCounts exec to config
#TODO: add temp() wrappers to bam and bai files
rule tagReads_withDir:
	input:
		TAR_GTF = '{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.withDir.refFlat.gtf'
	output:
		OUT_BAM='{DATADIR}/{sample}/TAR/STARsolo/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam'
		# TAR_BULK_COUNTS='{DATADIR}/{sample}/TAR/pseudobulk_featureCounts.txt'
	params:
		IN_BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
	threads:
		CORES
	log:
		'{DATADIR}/{sample}/TAR/dropseq_tag.log'
		# '{DATADIR}/{sample}/TAR/pseudobulk_featureCounts.txt.summary'
	shell:
		"""
		{DROPSEQ}/TagReadWithGeneFunction\
			I={params.IN_BAM}\
			O={output.OUT_BAM}\
			ANNOTATIONS_FILE={input.TAR_GTF}\
			GENE_NAME_TAG=XT\
			GENE_STRAND_TAG=GS\
			USE_STRAND_INFO=true\
			MAX_RECORDS_IN_RAM=10000000\
			CREATE_INDEX=true 2> {log}
		"""

rule sort_index_tagged_bam:
	input:
		IN_BAM = '{DATADIR}/{sample}/TAR/STARsolo/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam'
	output:
		OUT_BAM='{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bam',
		bamIndex='{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bam.bai'
	threads:
		CORES
	shell:
		"""
		samtools sort -@ {threads} {input.IN_BAM} -o {output.OUT_BAM}
		samtools index -@ {threads} {output.OUT_BAM}
		"""

# Get counts matrix for HMM-annotated features
rule extract_HMM_expression_withDir:
	input:
		BAM = '{DATADIR}/{sample}/TAR/TAR_tagged_withDir.bam'
	output:
		COUNTS='{DATADIR}/{sample}/TAR/TAR_expression_matrix_withDir.tsv.gz'
	params:
		BARCODES = '{DATADIR}/{sample}/STARsolo/Solo.out/GeneFull/filtered/barcodes.tsv.gz',
		TMPDIR = TMPDIR
	log:
		'{DATADIR}/{sample}/TAR/DropSeq_DigitalExpression.log'
	shell:
		"""
		{DROPSEQ}/DigitalExpression\
			I={input.BAM}\
			O={output.COUNTS}\
			TMP_DIR={params.TMPDIR}\
			CELL_BARCODE_TAG=CB\
			MOLECULAR_BARCODE_TAG=UB\
			GENE_NAME_TAG=XT\
			GENE_STRAND_TAG=GS\
			SUMMARY={log}\
			USE_STRAND_INFO=True\
			CELL_BC_FILE={params.BARCODES}\
			OMIT_MISSING_CELLS=False
		"""

# Might switch for umi-tools in the future
		# """
		# umi_tools count \
		# --per-gene \
		# --extract-umi-method=tag \
		# --assigned-status-tag=XS \
		# --gene-tag=XT \
		# --cell-tag=CB \
		# --umi-tag=UB \
		# --per-cell \
		# --wide-format-cell-counts \
		# --log={log} \
		# -I {input.BAM} \
		# -E {log} \
		# -S {output.COUNTS}
		# """

#########################################################
# Diff. TAR expression & BLAST
#########################################################

# generate differentially expressed genes and uTARs
rule getDiffFeatures:
	input:
		hmmFile='{DATADIR}/{sample}/TAR/TAR_expression_matrix_withDir.tsv.gz'
	output:
		diffFeatures='{DATADIR}/{sample}/TAR/diff_Expressed_Features.txt',
		MTX_DIR=directory('{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix'),
		MTX_BC='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/barcodes.tsv.gz',
		MTX_FEAT='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/features.tsv.gz',
		MTX_MTX='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/matrix.mtx.gz'
	params:
		geneFile='{DATADIR}/{sample}/STARsolo/Solo.out/GeneFull/filtered'
	shell:
		"""
		Rscript scripts/analyzeExpressionMat.R {params.geneFile} {input.hmmFile}
		"""

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

# extract diff uTAR fastas
rule getDiffSeqsToBlastFa:
	input:
		fastaRegions='{DATADIR}/{sample}/TAR/TAR_seqsToBlast.txt'
	output:
		TARfasta='{DATADIR}/{sample}/TAR/TAR_seqsToBlast.fa'
	params:
		fastaFile=GENOME_FASTA,
		fastaIndex=GENOME_FASTA+'.dict'
	shell:
		"""
		samtools faidx -r {input.fastaRegions} {params.fastaFile} > {output}
		"""

# run blast on fastas
rule ruleBlast:
	input:
		TARfasta='{DATADIR}/{sample}/TAR/TAR_seqsToBlast.fa'
	output:
		'{DATADIR}/{sample}/TAR/TAR_blastResults.txt'
	threads: 
		CORES
	params:
		blastDB=config['BLASTDB']
	shell:
		"""
		blastn -db {params.blastDB}/nt -query {input.TARfasta} -out {output} -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 5 -num_threads {CORES}
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
