# convert GTF to REFFlat, save in your cellranger reference
# 	not run if REFFlat file exists -  will only need to run this once for each reference genome
rule convertToRefFlat:
	input:
		GENES = GENES_GTF
	output:
		REFFLAT = GENES_GTF + ".refFlat"
	shell:
		"""
		{GTFTOGENEPRED} -genePredExt -geneNameAsName2 {input.GENES} refFlat.tmp
		paste <(cut -f 12 refFlat.tmp) <(cut -f 1-10 refFlat.tmp) > {output.REFFLAT}
		rm refFlat.tmp
		"""

# Run HMM
rule calcHMMbed:
	input:
		BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam' #need absolute path
	output:
		BED = temp('{DATADIR}/{sample}/TAR/TAR_reads.bed.gz')
	threads:
		CORES
	params:
		MEM = "64G" # Change this if you are using a low-RAM machine
	run:
		shell(
			f"""
			bash scripts/HMM.bash {input.BAM} {threads} {params.MEM} {MERGEBP} {THRESH} {wildcards.DATADIR}/{wildcards.sample}/TAR
			"""
		)

# Generate annotations from HMM in refFlat format
#		Returns annotations with and without considering strandedness/direction
rule calcHMMrefFlat:
	input:
		BEDGZ = '{DATADIR}/{sample}/TAR/TAR_reads.bed.gz',
		REFFLAT= GENES_GTF + ".refFlat"
	output:
		WITHDIR ='{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.refFlat'
	threads:
		CORES
	shell:
		"""
		Rscript scripts/generate_refFlat_script_both.R {input.REFFLAT} {input.BEDGZ} {threads}
		"""

# Convert stranded annotations to gtf format
rule HMM_refFlat_to_gtf:
	input:
		REFFLAT = '{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.refFlat'
	output:
		GTF = '{DATADIR}/{sample}/TAR/TAR_reads.bed.gz.refFlat.gtf'
	shell:
		"""
		Rscript scripts/convertRefFlatToGTF.R {input.REFFLAT}
		"""