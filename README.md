# uTAR_lite
Lightweight pipeline for discovering new features in single-cell and spatial transcriptomics data.

## Description
This is a lightweight version of the workflow for [TAR-scRNA-seq](https://github.com/fw262/TAR-scRNA-seq) that just generates a count matrix for unannotated transcriptionally active regions (uTARs). We removed a few of the dependencies and sped up the runtimes.

Right now, it only works for datasets which have already been aligned with [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md). All output files are saved into the output directory (```/path/to/sampleID/TAR```). 

## Example analyses:
#TODO
See `./vignettes` for example notebooks on how to use uTARs with your data 

## TODO:
- Add in .bam deduplication prior to HMM (saves run time, improves scalability, probably also improves sensitivity)
- Rewrite parameterization of THRESH (minnumber of UMIs, not relative to total reads)
- 

## Required Software
This workflow requires the following packages listed below. Please ensure that tool can be called from the command line (i.e. the paths to each tool is in your path variable). We recommend using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to organize your packages. We also included a .yml (```scTAR_STARsolo.yml```) file which can be used to initialize a conda environment with all of the required dependencies.

### 1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

### 2. [umi_tools](https://umi-tools.readthedocs.io/en/latest/index.html)
We replaced `DropSeqUtils` & `Picard` with `umi_tools`.

```
conda install -c bioconda -c conda-forge umi_tools==1.1.2
```

### 3. [R, version 3.6 or greater](https://www.r-project.org/)

Please also ensure that you have downloaded the following R packages. They will be used throughout the pipeline.
- [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [groHMM](https://www.bioconductor.org/packages/release/bioc/html/groHMM.html)
- [Seurat, version >= 4.0](https://satijalab.org/seurat/install.html)
- [data.table](https://github.com/Rdatatable/data.table)
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html)

### 4. [Samtools](http://www.htslib.org/)
```
conda install -c bioconda samtools
```

### 5. [GtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)
This tool is used to convert gtf annotation files to refFlat format.
```
conda install -c bioconda ucsc-gtftogenepred
```

### 6. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
Please make sure this tool is available in your working environment.
```
conda install -c bioconda bedtools
```

## Procedure

### 1. Align samples with STARsolo.
- Be sure that the aligned /bam files hav the `CB` and

### 2. Clone this repository.

Run the following command in your command line.
```
git clone git@github.com:mckellardw/uTAR_lite.git
```

### 3. Download required software listed above.

We recommend using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to manage packages, see above.

### 4. Align .fastq files with [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

Place each STAR output into the same parent directory, so that it follows this file tree: (`DATADIR` in `STARsolo_config.yaml`).
```
DATADIR
├── sample_A
│   ├── STARsolo
│   │   ├── Aligned.sortedByCoord.dedup.out.bam
│   │   ├── Aligned.sortedByCoord.dedup.out.bam.bai
│   │   ├── Aligned.sortedByCoord.out.bam
│   │   ├── Aligned.sortedByCoord.out.bam.bai
│   │   ├── Solo.out
│   │   │   ├── Gene
│   │   │   │   ├── filtered
│   │   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   │   ├── features.tsv.gz
│   │   │   │   │   ├── matrix.mtx.gz
├── sample_B
│   ├── STARsolo
│   │   ├── Aligned.sortedByCoord.dedup.out.bam
│   │   ├── Aligned.sortedByCoord.dedup.out.bam.bai
│   │   ├── Aligned.sortedByCoord.out.bam
│   │   ├── Aligned.sortedByCoord.out.bam.bai
│   │   ├── Solo.out
│   │   │   ├── Gene
│   │   │   │   ├── filtered
│   │   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   │   ├── features.tsv.gz
│   │   │   │   │   ├── matrix.mtx.gz
```

### 5. Edit the config.yaml file for your experiment.

Please change the variable names in the `config.yaml` as required for your analysis. This includes the following changes:
- **Samples**: STARsolo output directory names; each should be a directory in `DATADIR`
- **GENOME_FASTA**: path to the .fasta (genome sequence) file used to build STAR reference
- **GENOME_GTF**: path to the .gtf (gene annotations) file used to build STAR reference
- **DATADIR**: Path to where the STARsolo outputs (`.../DATADIR/sample_A`) are stored. Also note that outputs will be stored in each sample's directory, individually. See above for the expected file tree.
- **TMPDIR**: Directory to store temporary files.
- **GTFTOGENEPRED**: Path to gtfToGenePred tool. **NOTE** This step will only need to run once for each reference genome you use, and the resulting .refFlat file ("genes.refFlat") will be stored in the same directory as `GENOME_GTF`.
- **CORES**: Number of cores used in each step of the pipeline. To run multiple samples in parallel, please specify total number of cores in the snakemake command (i.e. "snakemake -j {total cores}").
- **MERGEBP**: Number of bases to merge in groHMM. Smaller numbers creates more TARs but takes longer to run. We recommend keeping the default value of 500.
- **THRESH**: Used to set TARs coverage threshold. This is sequence depth dependent. Default coverage threshold set at 1 in 10,000,000 uniquely aligned reads. For example, if there are 500,000,000 total aligned reads, TARs with at least 50 reads are kept when **THRESH** is set to 10,000,000. A higher **THRESH** value increases the number of TARs kept after filtering. We recommend keeping the default value of 10000000.

### 6. Run snakemake with the command ```snakemake -j [# total cores]```.

From the command line, `cd` into the directory ```.../TAR-SCRNA-seq/from_STARsolo```

Please ensure the Snakefile and STARsolo_config.yaml files as well as the scripts folder are in the directory where you intend to run the pipeline.

## Output files
#TODO- update output tree
All output files are stored in a directory inside of the STARsolo output  ```.../DATADIR/sample_A/TAR/```
- RefFlat format of TAR features with and without consideration of directionality stored in "**TAR_reads.bed.gz.withDir.genes.refFlat**" and "**TAR_reads.bed.gz.noDir.refFlat.refFlat**".
- Digital expression matrix for TAR features, with consideration of TAR directionality relative to annotated gene features, is stored in "**TAR_expression_matrix_withDir.txt.gz**". Additional rules within the **Snakefile** are provided to generate a matrix that does not consider strandedness- just uncomment them.
- A list of differentially expressed genes and uTARs in "**results_out/{sample}/{sample}\_diffMarkers.txt**".
- A list of differentially expressed uTARs and their labels based on BLASTn results in "**results_out/{sample}/{sample}\_TAR_diff_uTAR_Features_Labeled.txt**".
- Results of the BLASTn analysis for differentially expressed uTARs in "**TAR_blastResults.txt**".

