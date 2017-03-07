__author__ = 'alipirani'

import sys
import os
import argparse
import errno
from datetime import datetime
import ConfigParser
from config_settings import ConfigSectionMap
from modules.check_subroutines import *
from modules.stages import *
from modules.bedtools import *
from modules.qualimap import *
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp
from modules.logging_subprocess import *
from modules.log_modules import *
from modules.generate_PTR_dataframe import *

# Command Line Argument Parsing
def parser():
    parser = argparse.ArgumentParser(description='PTR Analysis pipeline for Illumina datasets.')
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-type', action='store', dest="type", help='Type of analysis: SE or PE', required=True)
    optional.add_argument('-config', action='store', dest="config", help='Path to Config file', required=False)
    required.add_argument('-PE1', action='store', dest="forward_raw", help='Path to Paired End file 1', required=True)
    optional.add_argument('-PE2', action='store', dest="reverse_raw", help='Path to Paired End file 2', required=False)
    required.add_argument('-o', action='store', dest="output_folder", help='Output Path ending with output directory name to save the results', required=True)
    required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique analysis name to save the results', required=True)
    required.add_argument('-index', action='store', dest="index", help='Reference Index Name. Change this argument in config file and mention the reference header name such as KP_NTUH_chr/KPNIH1/KPNIH32.', required=True)
    optional.add_argument('-c', action='store', dest="croplength", help='Crop Length in case needed')
    optional.add_argument('-f', action='store', dest="bam_input", help='Input Bam')
    return parser


# Main Pipeline
def pipeline(args, logger):

    keep_logging('\nSTART: Pipeline\n', 'START: Pipeline\n', logger, 'info')

    # Check Subroutines and create logger object: Arguments, Input files, Reference Index
    keep_logging('START: Checking Dependencies...', 'Checking Dependencies', logger, 'info')

    # Reference Genome file name
    reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
    keep_logging('Getting Reference Genome name from config file: {}'.format(reference), 'Getting Reference Genome name from config file: {}'.format(reference), logger, 'info')

    # Check FASTQ files
    if args.type != "PE":
        reverse_raw = "None"
        file_exists(args.forward_raw, args.forward_raw, reference)
    else:
        file_exists(args.forward_raw, args.reverse_raw, reference)

    # Check Java Version
    java_check()
    keep_logging('END: Checking Dependencies...\n', 'END: Checking Dependencies\n', logger, 'info')

    ## 1. Pre-Processing Raw reads using Trimmomatic
    keep_logging('START: Pre-Processing Raw reads using Trimmomatic', 'START: Pre-Processing Raw reads using Trimmomatic', logger, 'info')
    if args.type == "PE":
        trimmomatic(args.forward_raw, args.reverse_raw, args.output_folder, args.croplength, logger, Config)
    else:
        reverse_raw = "None"
        #trimmomatic(args.forward_raw, reverse_raw, args.output_folder, args.croplength, logger, Config)
    keep_logging('END: Pre-Processing Raw reads using Trimmomatic\n', 'END: Pre-Processing Raw reads using Trimmomatic\n', logger, 'info')

    # ## 2. Stages: Alignment using BWA
    # keep_logging('START: Mapping Reads using {}'.format(ConfigSectionMap("pipeline", Config)['aligner']), 'START: Mapping Reads using {}'.format(ConfigSectionMap("pipeline", Config)['aligner']), logger, 'info')
    # split_field = prepare_readgroup(args.forward_raw, ConfigSectionMap("pipeline", Config)['aligner'], logger)
    # files_to_delete = []
    # out_sam = align(args.bam_input, args.output_folder, args.index, split_field, args.analysis_name, files_to_delete, logger, Config, args.type)
    # keep_logging('END: Mapping Reads using {}\n'.format(ConfigSectionMap("pipeline", Config)['aligner']), 'END: Mapping Reads using {}\n'.format(ConfigSectionMap("pipeline", Config)['aligner']), logger, 'info')
    #
    # ## 3. Stages: Post-Alignment using SAMTOOLS, PICARD etc
    # keep_logging('START: Post-Alignment using SAMTOOLS, PICARD etc...', 'START: Post-Alignment using SAMTOOLS, PICARD etc...', logger, 'info')
    # out_sorted_bam = prepare_bam(out_sam, args.output_folder, args.analysis_name, files_to_delete, logger, Config)
    # # out_sorted_bam = "%s/%s_aln_sort.bam" % (args.output_folder, args.analysis_name)
    # final_coverage_file = bedtools(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
    # keep_logging('END: Post-Alignment using SAMTOOLS, PICARD etc...\n', 'END: Post-Alignment using SAMTOOLS, PICARD etc...\n', logger, 'info')
    #
    # ## 4. Stages: Statistics
    # keep_logging('START: Generating Statistics Reports', 'START: Generating Statistics Reports', logger, 'info')
    # alignment_stats_file = alignment_stats(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
    # gatk_DepthOfCoverage(out_sorted_bam, args.output_folder, args.analysis_name, reference, logger, Config)
    # qualimap_report = qualimap(out_sorted_bam, args.output_folder, args.analysis_name, logger, Config)
    # # final_coverage_file = "%s/%s_coverage.bed" % (args.output_folder, args.analysis_name)
    # keep_logging('END: Generating Statistics Reports\n', 'END: Generating Statistics Reports\n', logger, 'info')
    final_coverage_file = "%s/%s_coverage.bed" % (args.output_folder, args.analysis_name)
    ## 5. Stages: PTR Analysis
    keep_logging('START: Analyzing Bedfiles for PTR analysis', 'START: Analyzing Bedfiles for PTR analysis', logger, 'info')
    generate_PTR_dataframe(final_coverage_file, args.output_folder, logger, Config)
    keep_logging('END: Analyzing Bedfiles for PTR analysis\n', 'END: Analyzing Bedfiles for PTR analysis\n', logger, 'info')

## Check Subroutines
def usage():
    print "Usage: python pipeline.py [-h] -PE1 path-to-forward-PE-read -PE2 path-to-reverse-PE-read -o path-to-OUTPUT_FOLDER -analysis ANALYSIS_NAME -index INDEX_NAME_as_per_config_file \n"

# Validate Filenames for any unsupported characters
def Validate_filename( name ):
    pattern_strings = ['\.', '\&', '\>', 'aaa', '\*']
    pattern_string = '|'.join(pattern_strings)
    searchobj = re.search(pattern_string, name, flags=0)
    if searchobj:
        print "The file " + name + " contains unsupported characters such as quotes, spaces, or &:%?*><\$. \nPlease Provide another file name.\n"
        exit()

def file_exists(path1, path2, reference):

    if not os.path.isfile(path1):
        file_basename = os.path.basename(path1)
        keep_logging('The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()
    if path2 is not None:
        if not os.path.isfile(path2):
            file_basename = os.path.basename(path2)
            keep_logging('The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
            exit()
    if not os.path.isfile(reference):
        file_basename = os.path.basename(reference)
        keep_logging('The reference fasta file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), 'The reference fasta file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()
    if ConfigSectionMap("pipeline", Config)['aligner'] == "bwa":
        ref_index_suffix1 = reference + ".bwt"
        ref_index_suffix2 = reference + ".amb"
        ref_index_suffix3 = reference + ".ann"
        ref_index_suffix4 = reference + ".sa"
        ref_index_suffix5 = reference + ".pac"
    elif ConfigSectionMap("pipeline", Config)['aligner'] == "bowtie":
        ref_index_suffix1 = reference + ".1.bt2"
        ref_index_suffix2 = reference + ".2.bt2"
        ref_index_suffix3 = reference + ".3.bt2"
        ref_index_suffix4 = reference + ".4.ebwt"
        ref_index_suffix5 = reference + ".rev.1.bt2"
        ref_index_suffix6 = reference + ".rev.2.bt2"
    else:
        print "Please change the aligner section in config file."

        print "Different Aligner in config file"

    if not os.path.isfile(ref_index_suffix1):
        # file_basename = os.path.basename(reference)
        keep_logging('The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5), 'The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5), logger, 'warning')
        create_index(reference, ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5)
    else:
        keep_logging('Index file already exists.', 'Index file already exists.', logger, 'info')

    ############################################
    ref_fai_index = reference + ".fai"
    if not os.path.isfile(ref_fai_index):
        # file_basename = os.path.basename(reference)
        keep_logging('The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), 'The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), logger, 'warning')
        create_fai_index(reference, ref_fai_index)
    else:
        keep_logging('Samtools fai Index file already exists.', 'Samtools fai Index file already exists.', logger, 'info')
    ############################################
    dict_name = os.path.splitext(os.path.basename(reference))[0] + ".dict"
    if not os.path.isfile(ConfigSectionMap(args.index, Config)['ref_path'] + "/" + dict_name):
        keep_logging('The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name), 'The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name), logger, 'warning')
        picard_seqdict(dict_name, reference)
    else:
        keep_logging('The reference seq dict file required for GATK and PICARD exists.', 'The reference seq dict file required for GATK and PICARD exists.', logger, 'info')



def java_check():
    keep_logging('Checking Java Availability...', 'Checking Java Availability...', logger, 'info')
    jd = sp.check_output(["java", "-version"], stderr=sp.STDOUT)
    jd_version = jd.split('\n', 1)[0]
    if len(jd) < 1:
        keep_logging('Unable to find a java runtime environment. The pipeline requires java 6 or later.', 'Unable to find a java runtime environment. The pipeline requires java 6 or later.', logger, 'exception')
    else:
        keep_logging('Java Availability Check completed ...{}'.format(jd_version), 'Java Availability Check completed ...{}'.format(jd_version), logger, 'info')

def fileformat(file1, file2, final_out):
    print "Checking File format....\n"
    if not file1.endswith('.fastq.gz'):
        base = os.path.basename(file1)
        os.path.splitext(base)
        file_1 = os.path.splitext(base)[0]
        cmdstring = "gzip -d " + file1 + " > " + final_out + file_1
        print "Compressing input file " + base
        os.system(cmdstring)


    if not file2.endswith('.fastq.gz'):
        base = os.path.basename(file2)
        os.path.splitext(base)
        file_2 = os.path.splitext(base)[0]
        cmdstring = "gzip -d " + file2 + " > " + final_out + file_2
        print "Compressing input file " + base
        os.system(cmdstring)

def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

def create_index(reference,ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5):
    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    keep_logging('Creating Index of reference fasta file for {} aligner.'.format(aligner), 'Creating Index of reference fasta file for {} aligner'.format(aligner), logger, 'info')
    if aligner == "bwa":
        cmd = "%s %s %s" % (ConfigSectionMap("bwa", Config)['base_cmd'], ConfigSectionMap("bwa", Config)['index'], reference)
        keep_logging(cmd, cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')
    else:
        print "Different Aligner in config file"

def create_fai_index(reference, ref_fai_index):
    keep_logging('Creating FAI Index using Samtools.', 'Creating FAI Index using Samtools.', logger, 'info')
    cmd = "%s %s %s" % (ConfigSectionMap("samtools", Config)['base_cmd'], ConfigSectionMap("samtools", Config)['faiindex'], reference)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Samtools FAI Indexing step. Exiting.', 'Error in Samtools FAI Indexing step. Exiting.', logger, 'exception')
        sys.exit(1)


    if not os.path.isfile(ref_fai_index):
        keep_logging('The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), 'The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), logger, 'exception')
    else:
        keep_logging('Samtools Fai Index file created.', 'Samtools Fai Index file created.', logger, 'info')

def picard_seqdict(dict_name, reference):
    #dict_name = os.path.splitext(os.path.basename(reference_filename))[0] + ".dict"
    keep_logging('Creating Sequence Dictionary using Picard.', 'Creating Sequence Dictionary using Picard.', logger, 'info')
    cmd = "java -jar %s/%s/%s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s/%s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("picard", Config)['picard_bin'], ConfigSectionMap("picard", Config)['base_cmd'], reference, ConfigSectionMap(args.index, Config)['ref_path'], dict_name)
    keep_logging(cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Picard Sequence Dictionary creation step. Exiting.', 'Error in Picard Sequence Dictionary creation step. Exiting.', logger, 'exception')
        sys.exit(1)



###

# Main Method
if __name__ == '__main__':
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    args = parser().parse_args()
    global config_file
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    global logger
    if args.output_folder != '':
        args.output_folder += '/'
    make_sure_path_exists(args.output_folder)
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    pipeline(args, logger)
    keep_logging('End: Pipeline', 'End: Pipeline', logger, 'info')









# extract_mapped_reads = "/home/apirani/bin/samtools-1.2/samtools view -b -F 4 %s > %s_mapped.bam" % (out_sorted_bam, out_sorted_bam)
# extract_fastq = "/home/apirani/bin/bedtools2-master/bin/bedtools bamtofastq -i %s_mapped.bam -fq %s/%s_mapped.fastq" % (out_sorted_bam, args.output_folder, args.analysis_name)
# shuff_fastq = "paste <(cat %s/%s_mapped.fastq) | paste - - - - | shuf | awk -F\'\\t\' \'{OFS=\"\\n\"; print $1,$2,$3,$4 > \"%s/%s_mapped_shuff.fastq\"}\'" % (out_path, args.analysis_name, out_path, args.analysis_name)
# print shuff_fastq
# #gzip_fastq = "gzip %s/%s_mapped.fastq" % (args.output_folder, args.analysis_name)
# print extract_mapped_reads + "\n" + extract_fastq + "\n" + shuff_fastq
# os.system(extract_mapped_reads)
# os.system(extract_fastq)
# #os.system(shuff_fastq)