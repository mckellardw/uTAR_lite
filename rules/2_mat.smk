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
	threads:
		config['CORES']
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

rule 2_sort_index_tagged_bam:
	input:
		BAM = '{DATADIR}/{sample}/TAR/Aligned.sortedByCoord.dedup.out.bam.featureCounts.bam'
	output:
		BAM='{DATADIR}/{sample}/TAR/TAR_tagged.bam'
	threads:
		config['CORES']
	shell:
		"""
		samtools sort -@ {threads} {input.BAM} -o {output.BAM}
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
	threads:
		config['CORES']
	shell:
		"""
		pigz -p{threads} {input.COUNT_MTX}
		"""

