
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
		bash scripts/HMM.bash {input} {threads} {params.MEM} {MERGEBP} {THRESH} {wildcards.DATADIR}/{wildcards.sample}/TAR
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