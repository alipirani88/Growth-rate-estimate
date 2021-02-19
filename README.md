# Growth-rate-analysis

#### This pipeline calculates peak to through ratio (PTR; ratio of copy numbers at origin of replication to terminus) from  mapped reads coverage.

The algorithm follows the procedure as described in the publication: Growth dynamics of gut microbiota in health and disease inferred from single metagenomic samples http://science.sciencemag.org/content/349/6252/1101.long with a few minor changes.

#### The pipeline runs sequentially as follows:
***

> 1. Pre-Processing Raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
> 2. Read Alignment using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
> 3. Post-Alignment steps using [SAMTOOLS](http://samtools.sourceforge.net/), [GATK](https://software.broadinstitute.org/gatk/), [PICARD](https://broadinstitute.github.io/picard/)
> 4. Coverage graph analysis using [Bedtools](http://bedtools.readthedocs.io/en/latest/)
> 5. Binning the sequencing coverage and calculating PTR using In-house scripts

**Output**:
***

The pipeline generates various alignment and bed output files from different tools at different steps that can be used for manual inspection. The final PTR results can be found in:

- ***prefix_PTR.txt*** Estimated PTR value for the coverage graph
- ***prefix_perc_coverage_graph.R*** script to plot sequence coverage graph with PTR values

## Installation

The pipeline can be set up in two easy steps:

> 1. Clone the github directory onto your system.

```
git clone https://github.com/alipirani88/Growth-rate-analysis.git

```

> 2. Use Growth-rate-analysis_env.yml file to create conda environment.

Create two new environments - varcall and varcall_gubbins
```
conda env create -f Growth-rate-analysis/Growth-rate-analysis_env.yml -n Growth-rate-analysis

```

Check installation

```
conda activate Growth-rate-analysis

python Growth-rate-analysis/growth_rate.py -h
```
