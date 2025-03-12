# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
rule umitools_dedupBAM:
    input:
        WHITELIST="{DATADIR}/{sample}/bc/whitelist.txt",
        SORTEDBAM="{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam",
    output:
        DEDUPBAM="{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam",
    resources:
        threads=config["CORES"],
        mem_mb=8192
    log:
        log="{OUTDIR}/{sample}/dedup.log",
        err="{OUTDIR}/{sample}/dedup.err",
    shell:
        """
        bash scripts/split_dedup.sh {input.SORTEDBAM} {input.WHITELIST} {resources.threads} {output.DEDUPBAM} {DATADIR}/{wildcards.sample}/tmp/dedup \
        1> {log.log} \
        2> {log.err}
        """


# Index the deduplicated .bam file
rule umitools_indexDedupBAM:
    input:
        SORTEDBAM="{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam",
    output:
        BAI="{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai",
    resources:
        threads=config["CORES"],
        mem_mb=4096
    shell:
        """
        samtools index -@ {resources.threads} {input.SORTEDBAM}
        """
