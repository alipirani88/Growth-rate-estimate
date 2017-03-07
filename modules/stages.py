__author__ = 'alipirani'

from config_settings import ConfigSectionMap
from modules.bwa import align_bwa
from modules.samtools import *
from modules.picard import *
from modules.gatk import *
from modules.check_subroutines import *
import os
import gzip
import re
from modules.bowtie import *
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from modules.trim import *
from modules.qualimap import *

## Prepare ReadGroup option for BWA alignment
def prepare_readgroup(forward_read, aligner, logger):
    keep_logging('Preparing ReadGroup Info', 'Preparing ReadGroup Info', logger, 'info')
    samplename = os.path.basename(forward_read)
    if forward_read.endswith(".gz"):
        ###
        output = gzip.open(forward_read, 'rb')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        id_name = split_field[1]
        id_name = id_name.strip()
        if aligner == "bowtie":
            split_field = "--rg-id %s --rg SM:%s --rg LB:1 --rg PL:Illumina" % (split_field[1], samplename)
        elif aligner == "bwa":
            split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        else:
            exit()
        return split_field

    elif forward_read.endswith(".fastq"):
        ###
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field
        print "\n################## End: ReadGroup Preparation ##################\n"


    elif forward_read.endswith(".fq"):
        ###
        output = open(forward_read, 'r')
        firstLine = output.readline()
        split_field = re.split(r":",firstLine)
        split_field = "\"" + "@RG" + "\\tID:" + split_field[1] + "\\tSM:" + samplename + "\\tLB:1\\tPL:Illumina" + "\""
        return split_field
        print "\n################## End: ReadGroup Preparation ##################\n"




## Raw data Pre-processing using Trimmomatic ##
def trimmomatic(input1, input2, out_path, crop, logger, Config):
    trim(input1, input2, out_path, crop, logger, Config)
## End ##

## bwa, smalt, bowtie: Alignment ##
def align(bam_input, out_path, ref_index, split_field, analysis, files_to_delete, logger, Config, type):
    reference = ConfigSectionMap(ref_index, Config)['ref_path'] + "/" + ConfigSectionMap(ref_index, Config)['ref_name']
    forward_clean = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['f_p']
    reverse_clean = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['r_p']
    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    if aligner == "bwa":
        base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bwa", Config)['bwa_bin'] + "/" + ConfigSectionMap("bwa", Config)['base_cmd']
        #check if the input is bam or fastq
        if bam_input:
            if bam_input.endswith(".bam"):
                #do alignment of bam and all here
                print "bam alignment"
            else:
                #throw error
                print "error"
        else:
            out_file = align_bwa(base_cmd,forward_clean, reverse_clean, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type)
            return out_file
    elif aligner == "smalt":
        print "Smalt addition pending"
        exit()
        usage()
    elif aligner == "bowtie":
        base_cmd = ConfigSectionMap("bin_path")['binbase'] + "/" + ConfigSectionMap("bowtie")['bowtie_bin'] + "/" + ConfigSectionMap("bowtie")['align_cmd']
        if bam_input:
            if bam_input.endswith(".bam"):
                #do alignment of bam and all here
                print "bam alignment"
            else:
                #throw error
                print "error"
        else:
            out_sam = align_bowtie(base_cmd,forward_clean, reverse_clean, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type)
            return out_sam

## bwa, smalt, bowtie: Alignment ##
def align(bam_input, out_path, ref_index, split_field, analysis, files_to_delete, logger, Config, type):
    reference = ConfigSectionMap(ref_index, Config)['ref_path'] + "/" + ConfigSectionMap(ref_index, Config)['ref_name']
    forward_clean = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['f_p']
    reverse_clean = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['r_p']
    forward_unpaired = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['f_up']
    reverse_unpaired = out_path + "/" + ConfigSectionMap("Trimmomatic", Config)['r_up']
    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    if aligner == "bwa":
        base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bwa", Config)['bwa_bin'] + "/" + ConfigSectionMap("bwa", Config)['base_cmd']
        #check if the input is bam or fastq
        if bam_input:
            if bam_input.endswith(".bam"):
                #do alignment of bam and all here
                print "bam alignment"
            else:
                #throw error
                print "error"
        else:
            out_file = align_bwa(base_cmd,forward_clean, reverse_clean, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type)
            return out_file
    elif aligner == "smalt":
        print "Smalt addition pending"
        exit()
        usage()
    elif aligner == "bowtie":
        base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bowtie", Config)['bowtie_bin'] + "/" + ConfigSectionMap("bowtie", Config)['align_cmd']
        if bam_input:
            if bam_input.endswith(".bam"):
                #do alignment of bam and all here
                print "bam alignment"
            else:
                #throw error
                print "error"
        else:
            out_sam = align_bowtie(base_cmd,forward_clean, reverse_clean, forward_unpaired, reverse_unpaired, out_path, reference, split_field, analysis, files_to_delete, logger, Config, type)
            return out_sam
## End ##

## samtools: Post-Alignment SAM/BAM conversion, sort, index ##
def prepare_bam(out_sam, out_path, analysis, files_to_delete, logger, Config):
    out_bam = samtobam(out_sam, out_path, analysis, files_to_delete, logger, Config)
    out_sort_bam = sort_bam(out_bam, out_path, analysis, logger, Config)
    files_to_delete.append(out_sort_bam)
    out_marked_bam = markduplicates(out_sort_bam, out_path, analysis, files_to_delete, logger, Config)
    out_sort_bam = sort_bam(out_marked_bam, out_path, analysis, logger, Config)
    index_bam(out_sort_bam, out_path, logger, Config)
    if not os.path.isfile(out_sort_bam):
        keep_logging('Error in SAM/BAM conversion, sort, index. Exiting.', 'Error in SAM/BAM conversion, sort, index. Exiting.', logger, 'exception')
        exit()
    else:
        return out_sort_bam
## End ##


## Statistics Report ##
def alignment_stats(out_sorted_bam, out_path, analysis, logger, Config):
    alignment_stats_file = flagstat(out_sorted_bam, out_path, analysis, logger, Config)
    keep_logging('The Alignments Stats file from Samtools: {}'.format(alignment_stats_file), 'The Alignments Stats file from Samtools: {}'.format(alignment_stats_file), logger, 'debug')
    return alignment_stats_file

def qualimap(out_sorted_bam, out_path, analysis, logger, Config):
    qualimap_report = bamqc(out_sorted_bam, out_path, analysis, logger, Config)
    return qualimap_report
## END ##


## Remove SAM files ##
def remove_files(analysis, out_path, out_sam, out_sorted_bam):
    os.remove(out_sam)
    if os.path.isfile(out_sorted_bam):
        raw_bam_file = "%s/%s_aln.bam" % (out_path, analysis)
        os.remove(raw_bam_file)
## END ##

