# Growth-rate-analysis (In-progress)

#### This pipeline estimates the growth rate of a bacteria in a medium by looking at their sequencing read coverage and calculating their peak to through ratio(PTR; ratio of copy numbers at origin of replication to terminus). 

The algorithm of this pipeline is based on a published article named: Growth dynamics of gut microbiota in health and disease inferred from single metagenomic samples http://science.sciencemag.org/content/349/6252/1101.long with a few minor changes.

#### The pipeline runs sequentially as follows:
***

> 1. Pre-Processing Raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
> 2. Read Alignment using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
> 3. Post-Alignment steps using [SAMTOOLS](http://samtools.sourceforge.net/), [GATK](https://software.broadinstitute.org/gatk/), [PICARD](https://broadinstitute.github.io/picard/) etc
> 4. Coverage graph analysis using [Bedtools](http://bedtools.readthedocs.io/en/latest/)
> 5. Binning the sequencing coverage and calculating PTR using In-house scripts

**Output**:
***

The pipeline generates various output files from different tools at different steps. The most notable ones are:

- ***Clean reads***: *.fq.gz files from trimmomatic.

- ***Alignment files***: analysisname_aln.sam and analysisname_aln.bam from Bowtie2, analysisname_aln_marked.bam from GATK MarkDuplicates, and finally a sorted BAM from marked bam file analysisname_aln_sort.bam. 

- ***Bed file***: Bedtools genome coverage file analysisname_coverage.bed

- ***coverage_graph.R*** script to plot sequence coverage graph with PTR values

Example:

<<<<<<< HEAD

=======
![](/perc_coverage_graph.png?raw=true)
>>>>>>> 9e55eedb374d38fa24b7e41530494eb68eb8cd0f



