# convert GTF to REFFlat, save in your cellranger reference
# 	not run if REFFlat file exists -  will only need to run this once for each reference genome
rule 1_convertToRefFlat:
    input:
        GENES = GENES_GTF
    output:
        REFFLAT = GENES_GTF.replace(".gtf", ".refFlat")
    conda:
        f"{workflow.basedir}/envs/ucsc.yml"
    resources:
        threads=1,
        mem_mb=2048
    shell:
        """
        gtfToGenePred -genePredExt -geneNameAsName2 {input.GENES} refFlat.tmp
        paste <(cut -f 12 refFlat.tmp) <(cut -f 1-10 refFlat.tmp) > {output.REFFLAT}
        rm refFlat.tmp
        """

# Run HMM
rule 1_calcHMMbed:
    input:
        BAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam' #need absolute path
    output:
        BED = temp('{DATADIR}/{sample}/TAR/TAR_reads.bed.gz')
    resources:
        threads=config['CORES'],
        mem_mb=65536
    params:
        MEM = "64G" 
    conda:
        f"{workflow.basedir}/envs/hmm.yml"
    shell:
        """
        bash scripts/HMM.bash \
            {input.BAM} \
            {resources.threads} \
            {params.MEM} \
            {MERGEBP} \
            {THRESH} \
            {wildcards.DATADIR}/{wildcards.sample}/TAR
        """

# Generate annotations from HMM in refFlat format
#		Returns annotations with and without considering strandedness/direction
rule 1_calcHMMrefFlat:
    input:
        BEDGZ = '{DATADIR}/{sample}/TAR/TAR_reads.bed.gz',
        REFFLAT= GENES_GTF.replace(".gtf", ".refFlat")
    output:
        WITHDIR = '{DATADIR}/{sample}/TAR/TAR_reads.withDir.refFlat'
    resources:
        threads=config['CORES'],
        mem_mb=8192
    log:
        '{DATADIR}/{sample}/TAR/calcHMMrefFlat.log'
    conda:
        f"{workflow.basedir}/envs/hmm.yml"
    shell:
        """
        Rscript scripts/generate_refFlat_script_both.R {input.REFFLAT} {input.BEDGZ} {resources.threads} | tee {log}
        """

# Convert stranded annotations to gtf format
rule 1_HMM_refFlat_to_gtf:
    input:
        REFFLAT = '{DATADIR}/{sample}/TAR/TAR_reads.withDir.refFlat'
    output:
        GTF = '{DATADIR}/{sample}/TAR/TAR_reads.withDir.gtf'
    resources:
        threads=1,
        mem_mb=2048
    conda:
        f"{workflow.basedir}/envs/hmm.yml"
    shell:
        """
        Rscript scripts/convertRefFlatToGTF.R {input.REFFLAT}
        """