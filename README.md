# Growth-rate-analysis

##### Growth-rate-analysis is a python pipeline to calculate peak to through ratio (PTR - ratio of copy numbers at origin of replication to terminus) from  mapped reads coverage.

The pipeline follows the coverage smoothing algorithm as described in the publication: 

[Growth dynamics of gut microbiota in health and disease inferred from single metagenomic samples](http://science.sciencemag.org/content/349/6252/1101.long)

## Installation

The pipeline can be set up in two easy steps:

> 1. Clone the github directory onto your system.

```
git clone https://github.com/alipirani88/Growth-rate-analysis.git

```

> 2. Use Growth-rate-analysis_env.yml file provided with the github repo to create an environment - Growth-rate-analysis.


```
conda env create -f Growth-rate-analysis/Growth-rate-analysis_env.yml -n Growth-rate-analysis

```

Check installation

```
conda activate Growth-rate-analysis

python Growth-rate-analysis/growth_rate.py -h
```

## Usage

Lets say you want to estimate PTR from the sequencing reads CFT073_condition_R1.fastq.gz. The pipeline will map the reads against CFT073 reference genome and calculate PTR from its smoothed coverage. 

```
python -type SE -PE1 /Path-to-forward-end-fastq-file/CFT073_condition_R1.fastq.gz -o /path-to-output-directory/ -analysis CFT073_condition -index CFT073 -steps All

```

- The above command will run the pipeline on fastq reads provided with -PE1 argument
- The results will be saved in the output directory with a prefix CFT073_condition.
- The reference genome and its path will be detected from the CFT073 settings that is set in config file.

#### config

config file is a high level easy to write YAML format configuration file that lets you configure your system wide runs and specify analysis parameters, path to the installed tools, data and system wide information.

- This config file will contain High level information such as locations of installed programs, cores and memory usage for running on HPC compute cluster, path to a reference genome, various parameters used by different tools. These settings will apply across multiple runs and samples. 

- The config file stores data in KEY: VALUE pair. 

- An example [config](https://github.com/alipirani88/Growth-rate-analysis/blob/master/config) file with default parameters is included with the installation folder. You can customize this config file and provide it with the -config argument or edit this config file based on your requirements. 

- If you wish to run pipeline in hpc compute environment such as PBS or SLURM, change the number of nodes/cores memory reuirements based on your needs else the pipeline will run with default settings.

#### index

- Add reference genome index name and its path to config file. For example; if you have set the reference genome path in config file as shown below, then the required value for command line argument -index would be -index CFT073

```
#index name
[CFT073]
# path to the reference genome fasta file.
Ref_Path: /nfs/esnitkin/bin_group/variant_calling_bin/reference/CFT073/
# Name of reference genome fasta file.
Ref_Name: CFT073.fasta
```

#### The pipeline runs sequentially as follows:
***

> 1. Pre-Processing Raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
> 2. Read Alignment using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
> 3. Post-Alignment steps using [SAMTOOLS](http://samtools.sourceforge.net/) and [PICARD](https://broadinstitute.github.io/picard/) duplicate reads removal
> 4. Coverage graph analysis using [Bedtools](http://bedtools.readthedocs.io/en/latest/)
> 5. Binning the sequencing coverage and calculating PTR using In-house script.

**Output**:
***

The pipeline generates various alignment and bed output files from different tools at different steps that can be used for manual inspection. The final PTR results can be found in:

- ***prefix_PTR.txt*** Estimated PTR value for the coverage graph
- ***prefix_perc_coverage_graph.R*** script to plot sequence coverage graph with PTR values


