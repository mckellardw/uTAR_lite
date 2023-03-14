########################################################################################################
## Tagging and labeling
##	source: https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html
########################################################################################################


# Label .bam file with each HMM feature
#TODO: add featureCounts exec to config
#TODO: add temp() wrappers to bam and bai files
# rule tagReads:
# 	input:
# 		BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam',
# 		TAR_GTF = '{DATADIR}/{sample}/TAR/TAR_reads.withDir.gtf'
# 	output:
# 		BAM='{DATADIR}/{sample}/TAR/STARsolo/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam'
# 	# params:
# 	# 	IN_BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
# 	threads:
# 		CORES
# 	log:
# 		'{DATADIR}/{sample}/TAR/dropseq_tag.log'
# 	run:
# 		shell(
# 			f"""
# 			{DROPSEQ_EXEC}/TagReadWithGeneFunction\
# 				I={input.BAM}\
# 				O={output.BAM}\
# 				ANNOTATIONS_FILE={input.TAR_GTF}\
# 				GENE_NAME_TAG=XT\
# 				GENE_STRAND_TAG=GS\
# 				USE_STRAND_INFO=true\
# 				MAX_RECORDS_IN_RAM=10000000\
# 				CREATE_INDEX=true 2> {log}
# 			"""
# 		)
		#Old command which used DropSeqUtils
		# f"""
		# {DROPSEQ_EXEC}/TagReadWithGeneFunction\
		# 	I={input.BAM}\
		# 	O={output.BAM}\
		# 	ANNOTATIONS_FILE={input.TAR_GTF}\
		# 	GENE_NAME_TAG=XT\
		# 	GENE_STRAND_TAG=GS\
		# 	USE_STRAND_INFO=true\
		# 	MAX_RECORDS_IN_RAM=10000000\
		# 	CREATE_INDEX=true 2> {log}
		# """


rule tagReads:
	input:
		BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam',
		TAR_GTF = '{DATADIR}/{sample}/TAR/TAR_reads.withDir.gtf'
	output:
		BAM=temp('{DATADIR}/{sample}/TAR/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam')
	# params:
	# 	IN_BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
	threads:
		CORES
	log:
		'{DATADIR}/{sample}/TAR/dropseq_tag.log'
	run:
		shell(
			f"""
			featureCounts \
			-T {threads} \
			-t exon \
			-g gene_id \
			-a {input.TAR_GTF} \
			--largestOverlap \
			--readExtension5 0 \
			--readExtension3 0 \
			-s 1 \
			-M \
			-o {DATADIR}/{wildcards.sample}/TAR/gene_assigned \
			-R BAM \
			{input.BAM} \
			2> {log}
			"""
		)
			# -o {output.BAM} \

rule sort_index_tagged_bam:
	input:
		IN_BAM = '{DATADIR}/{sample}/TAR/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam'
	output:
		OUT_BAM='{DATADIR}/{sample}/TAR/TAR_tagged.bam',
		bamIndex='{DATADIR}/{sample}/TAR/TAR_tagged.bam.bai'
	threads:
		CORES
	shell:
		"""
		{SAMTOOLS_EXEC} sort -@ {threads} {input.IN_BAM} -o {output.OUT_BAM}
		{SAMTOOLS_EXEC} index -@ {threads} {output.OUT_BAM}
		"""

# Get counts matrix for HMM-annotated features
rule extract_HMM_expression:
	input:
		BAM = '{DATADIR}/{sample}/TAR/TAR_tagged.bam'
	output:
		# COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR_matrix.mtx'
		COUNT_MTX=temp('{DATADIR}/{sample}/TAR/counts.tsv.gz')
	params:
		BARCODES = '{DATADIR}/{sample}/STARsolo/Solo.out/GeneFull/filtered/barcodes.tsv.gz'
		# TMPDIR = TMPDIR
	log:
		# '{DATADIR}/{sample}/TAR/DropSeq_DigitalExpression.log'
		'{DATADIR}/{sample}/TAR/umitools_count.log'
	run:
		shell(
			f"""
			{UMITOOLS_EXEC} count \
			--per-gene \
			--extract-umi-method=tag \
			--assigned-status-tag=XS \
			--gene-tag=XT \
			--cell-tag=CB \
			--umi-tag=UB \
			--per-cell \
			--log={log} \
			--stdin {input.BAM} \
			--log {log} \
			--stdout {output.COUNT_MTX}
			"""
		)
			# --multimapping-detection-method=NH \
			# --wide-format-cell-counts \
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

# Convert the long-format matrix (from umi_tools) to an .mtx file
rule counts_long2mtx:
	input:
		COUNT_MTX='{DATADIR}/{sample}/TAR/counts.tsv.gz'
	output:
		COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR.mtx'
	run:
		shell(
			"""
			python scripts/long2mtx.py {input.COUNT_MTX} {output.COUNT_MTX} --output-format mtx
			"""
		)
		
# Compress the count matrix
rule gzip_counts:
	input:
		COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR.mtx'
	output:
		COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR.mtx.gz'
	threads:
		CORES
	shell:
		"""
		pigz -p{threads} {input.COUNT_MTX}
		"""

# Convert the long-format matrix (from umi_tools) to an .h5 file
rule counts_long2h5:
	input:
		COUNT_MTX='{DATADIR}/{sample}/TAR/counts.tsv.gz'
	output:
		COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR.h5'
	run:
		shell(
			"""
			python scripts/long2mtx.py {input.COUNT_MTX} {output.COUNT_MTX} --output-format h5
			"""
		)

# generate differentially expressed genes and uTARs
# rule getDiffFeatures:
# 	input:
# 		hmmFile='{DATADIR}/{sample}/TAR/TAR_expression_matrix.tsv.gz'
# 	output:
# 		diffFeatures='{DATADIR}/{sample}/TAR/diff_Expressed_Features.txt',
# 		MTX_DIR=directory('{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix'),
# 		MTX_BC='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/barcodes.tsv.gz',
# 		MTX_FEAT='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/features.tsv.gz',
# 		MTX_MTX='{DATADIR}/{sample}/TAR/TAR_feature_bc_matrix/matrix.mtx.gz'
# 	params:
# 		geneFile='{DATADIR}/{sample}/STARsolo/Solo.out/GeneFull/filtered'
# 	shell:
# 		"""
# 		Rscript scripts/analyzeExpressionMat.R {params.geneFile} {input.hmmFile}
# 		"""