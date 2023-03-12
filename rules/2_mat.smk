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
	params:
		IN_BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
	threads:
		CORES
	log:
		'{DATADIR}/{sample}/TAR/dropseq_tag.log'
	shell:
		"""
		{DROPSEQ_EXEC}/TagReadWithGeneFunction\
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
		# '{DATADIR}/{sample}/TAR/DropSeq_DigitalExpression.log'
        '{DATADIR}/{sample}/TAR/umitools_count.log'
	shell:
        """
		{UMITOOLS_EXEC} count \
		--per-gene \
		--extract-umi-method=tag \
		--assigned-status-tag=XS \
		--gene-tag=XT \
		--cell-tag=CB \
		--umi-tag=UB \
		--per-cell \
		--wide-format-cell-counts \
		--log={log} \
		-I {input.BAM} \
		-E {log} \
		-S {output.COUNTS}
		"""

		# """
		# {DROPSEQ}/DigitalExpression\
		# 	I={input.BAM}\
		# 	O={output.COUNTS}\
		# 	TMP_DIR={params.TMPDIR}\
		# 	CELL_BARCODE_TAG=CB\
		# 	MOLECULAR_BARCODE_TAG=UB\
		# 	GENE_NAME_TAG=XT\
		# 	GENE_STRAND_TAG=GS\
		# 	SUMMARY={log}\
		# 	USE_STRAND_INFO=True\
		# 	CELL_BC_FILE={params.BARCODES}\
		# 	OMIT_MISSING_CELLS=False
		# """
	


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