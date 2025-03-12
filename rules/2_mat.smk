########################################################################################################
## Tagging and labeling
##	source: https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html
########################################################################################################


# Label .bam file with each HMM feature
rule 2_tagReads:
	input:
		BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam',
		TAR_GTF = '{DATADIR}/{sample}/TAR/TAR_reads.withDir.gtf'
	output:
		BAM=temp('{DATADIR}/{sample}/TAR/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam')
	resources:
		threads=config['CORES'],
		mem_mb=8192
	log:
		'{DATADIR}/{sample}/TAR/dropseq_tag.log'
	shell:
		"""
		featureCounts \
			-T {resources.threads} \
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

rule 2_sort_index_tagged_bam:
	input:
		BAM = '{DATADIR}/{sample}/TAR/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam'
	output:
		BAM='{DATADIR}/{sample}/TAR/TAR_tagged.bam'
	resources:
		threads=config['CORES'],
		mem_mb=4096
	shell:
		"""
		samtools sort -@ {resources.threads} {input.BAM} -o {output.BAM}
		"""

# Get counts matrix for HMM-annotated features
rule 2_extract_HMM_expression:
	input:
		BAM = '{DATADIR}/{sample}/TAR/TAR_tagged.bam'
		BAI = '{DATADIR}/{sample}/TAR/TAR_tagged.bam.bai'
	output:
		# COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR_matrix.mtx'
		COUNT_MTX='{DATADIR}/{sample}/TAR/counts.tsv.gz'
	params:	
        CELL_TAG="CB",  # uncorrected = CR
        GENE_TAG="XT",  #GN XS
        UMI_TAG="UB",
		STATUS_TAG="XS"
	log:
		log='{DATADIR}/{sample}/TAR/umitools_count.log',
		err='{DATADIR}/{sample}/TAR/umitools_count.err'
    conda:
        f"{workflow.basedir}/envs/umi_tools.yml"
	resources:
		threads=config['CORES'],
		mem_mb=8192
	shell:
		"""
		umi_tools count --extract-umi-method=tag \
			--per-cell \
			--per-gene \
			--assigned-status-tag={params.CELL_TAG} \
			--cell-tag={params.CELL_TAG} \
			--gene-tag={params.GENE_TAG}  \
			--umi-tag={params.UMI_TAG}  \
			--log={log.log} \
			--stdin {input.BAM} \
			--stdout {output.COUNT_MTX} \
		2> {log.err}
		"""
			# --multimapping-detection-method=NH \

# Convert the long-format matrix (from umi_tools) to an .mtx file
rule 2_counts_long2mtx:
	input:
		TSV='{DATADIR}/{sample}/TAR/counts.tsv.gz'
	output:
		MTX='{DATADIR}/{sample}/TAR/uTAR.mtx',
		GENES='{DATADIR}/{sample}/TAR/uTAR_genes.tsv.gz',
		CELLS='{DATADIR}/{sample}/TAR/uTAR_cells.tsv.gz'
	resources:
		threads=1,
		mem_mb=2048
	shell:
		"""
		python scripts/long2mtx.py \
			--umitools_tsv {input.TSV} \
			--out_mat {output.MTX} \
			--output-format mtx
		"""
		
# Compress the count matrix
rule 2_gzip_counts:
	input:
		COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR.mtx'
	output:
		COUNT_MTX='{DATADIR}/{sample}/TAR/uTAR.mtx.gz'
	resources:
		threads=config['CORES'],
		mem_mb=2048
	shell:
		"""
		pigz -p{resources.threads} {input.COUNT_MTX}
		"""

