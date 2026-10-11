# NinaBee13/cw-rnaseq-roehre-schuele

## Introduction

**NinaBee13/cw-rnaseq-roehre-schuele** is a bioinformatics pipeline that ingests RNA sequencing data, and provides the user with a comprehensive summary of the processing steps within the pipeline. The input should have the form of a samplesheet, with FASTQ files, strandedness, and sample ids. This input is than processed by a series of rna-sep processing steps, such as Quality Control (FASTQC), trimming (FASTP), pseudo-/alignment (HISAT2, Salmon), marking of duplicate reads (PICARD).

<!-- TODO nf-core:
Thi
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->
The pipeline accepts single-end or paired-end FASTQ files, including multiple
sequencing runs per biological sample. Reads are merged per sample, quality
controlled and trimmed before being processed in two parallel analysis branches.
HISAT2 is used for genome alignment and Salmon for transcript/gene abundance
estimation. Alignment and quality-control statistics are summarized with MultiQC.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2, strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz, 
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->
### Samplesheet

The input samplesheet must contain the following columns:

```csv
sample,fastq_1,fastq_2,strandedness
WT_REP1,/path/run1_R1.fastq.gz,/path/run1_R2.fastq.gz,reverse
WT_REP1,/path/run2_R1.fastq.gz,/path/run2_R2.fastq.gz,reverse
WT_REP2,/path/WT_REP2_R1.fastq.gz,/path/WT_REP2_R2.fastq.gz,reverse

Now, you can run the pipeline using:

```bash
nextflow run NinaBee13/cw-rnaseq-roehre-schuele \
    -profile docker \
    --input samplesheet.csv \
    --fasta reference.fa \
    --gtf annotation.gtf \
    --outdir results

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

For a minimal example run:

```bash
nextflow run NinaBee13/cw-rnaseq-roehre-schuele \
    -profile test,docker \
    --outdir test_results

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

NinaBee13/cw-rnaseq-roehre-schuele was originally written by Nina Röhre, Laura Schüle.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use NinaBee13/cw-rnaseq-roehre-schuele for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
