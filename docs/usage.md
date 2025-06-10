# Zymo-research/aladdin-viralrecon (illumina): Usage

This document describes how to use [Aladdin Viralrecon (Illumina) pipeline](https://github.com/Zymo-Research/aladdin-viralrecon)

## Table of contents
1.  [Input and Output Options](#input-and-output-options)
2.  [Main options](#main-options)
3.  [FastQ Preprocessing](#fastq-preprocessing)
4.  [QC, read trimming and filtering options](#illumina-qc-read-trimming-and-filtering-options)
5.  [Variant calling options](#illumina-variant-calling-options)
6.  [De novo assembly](#illumina-de-novo-assembly-options)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in ~/.bashrc or ~./bash_profile):

```
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the full pipeline is as follows:

```
nextflow run Zymo-Research/aladdin-viralrecon \
  --project "<your project name>" \
  --genome "MN908947.3" \
  --design "<your design file>" \
  --viral_db "<your viral db>" \
  --kraken2_db "<your kraken2 db>" \
  -profile <docker/singularity/...> \
  -work-dir "<work dir>" \
  --outdir "<output dir>" \
  -r "main" 
```

This command will retrieve the pipeline code from GitHub. This requires that you have setup appropriate GitHub access privilidges. If you have downloaded the pipeline code, you can subsitute `Zymo-Research/aladdin-viralrecon` with `main.nf`.

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```
nextflow pull Zymo-Research/aladdin-viralrecon
```

## Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Zymo-Research/aladdin-viralrecon](https://github.com/Zymo-Research/aladdin-viralrecon) releases page and find the latest version number - numeric only (eg. 0.1.0). Then specify this when running the pipeline with -r (one hyphen) - eg. -r 0.1.0.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Input and Output Options:
- **--design**:
  You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Path to a comma-separated text file (CSV) containing information about the samples in the experiment. The CSV file must have a header row with the following columns:
  - sample
  - read_1
  - read_2

  It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.
  ```csv title="samplesheet.csv"
  sample,fastq_1,fastq_2
  SAMPLE_1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
  SAMPLE_1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
  SAMPLE_2,AEG588A2_S4_L003_R1_001.fastq.gz,
  ```

- **--outdir**:
  Directory where the results will be saved. Use absolute paths, especially if you are working on cloud infrastructure.

## Reference genome options:

- **--genome**:
  The name of the viral reference genome you want to align to. For more information regarding keys, please consult the [Genomes config file](https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config).


## Illumina QC, read trimming and filtering options:

- **--kraken2_variants_host_filter**: 
  You can choose to remove host reads identified by Kraken2 before running variant calling steps in the pipeline.

- **--kraken2_assembly_host_filter**:
  You can choose to remove host reads identified by Kraken2 before running de novo assembly steps in the pipeline.

## Illumina variant calling options:

- **--min_mapped_reads**:
  Minimum number of mapped reads below which samples are removed from further processing. Some downstream steps in the pipeline will fail if this threshold is too low.

- **--ivar_trim_noprimer**:
  You can choose this option to discard reads without primers for ivar trim.

- **--ivar_trim_offset**:
  By selecting an offset, tye reads that occur at the specified offset positions relative to primer positions will also be trimmed. Please note that this option will need to be set for some amplicon-based sequencing protocols (e.g. SWIFT) as described and implemented [here](https://github.com/andersen-lab/ivar/pull/88).

 ## Illumina de novo assembly options:

- **--spades_mpde**:
  You can specify the SPAdes mode you would like to run (default: `rnaviral`).

- **--min_contig_length**:
  For this option, the minim contig length will be filtered from BLAST results.

- **--min_perc_contig_aligned**:
  The minimum percentage of contig aligned to filter from BLAST results.
