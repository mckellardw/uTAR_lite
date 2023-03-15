# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
rule umitools_dedupBAM:
    input:
        WHITELIST = "{DATADIR}/{sample}/bb/whitelist.txt", # TODO update this...
        SORTEDBAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        DEDUPBAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam' 
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/dedup.log'
    run:
        shell(
            f"""
            bash scripts/split_dedup.sh {input.SORTEDBAM} {input.WHITELIST} {threads} {output.DEDUPBAM} {DATADIR}/{wildcards.sample}/tmp/dedup | tee {log}
            """
        )

# Index the deduplicated .bam file
rule umitools_indexDedupBAM:
    input:
        SORTEDBAM = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam'
    output:
        BAI = '{DATADIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai'
    threads:
        config['CORES']
    shell:
        """
        {SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
        """