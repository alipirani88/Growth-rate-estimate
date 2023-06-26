__author__ = 'alipirani'
import os
import sys
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp
import errno
import gzip
import re
import subprocess
from modules.log_modules import keep_logging
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *

def trim(input1, input2, out_path, crop, adapter, logger, Config):
    if str(adapter) == "None":
        proc = subprocess.Popen(["which trimmomatic"], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        adapter_file = (out.strip()).replace('bin', 'share') + "/adapters/NexteraPE-PE.fa"
    else:
        adapter_file = adapter

    if input2 != "None":
        keep_logging('Pre-processing PE reads using Trimmomatic.', 'Pre-processing PE reads using Trimmomatic.', logger, 'info')
        #adapter_file = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "/" + ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath']
        clean_filenames = out_path + ConfigSectionMap("Trimmomatic", Config)['f_p'] + " " + out_path + ConfigSectionMap("Trimmomatic", Config)['f_up'] + " " + out_path + ConfigSectionMap("Trimmomatic", Config)['r_p'] + " " + out_path + ConfigSectionMap("Trimmomatic", Config)['r_up']
        # changing this parameter for KPC variant analysis for keeping both reads. date: 31 August
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] +  ConfigSectionMap("Trimmomatic", Config)['minadapterlength'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['keep_both_reads']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']
        if not crop:
            # cmdstring = "java -jar " + ConfigSectionMap("bin_path", Config)['binbase'] + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "trimmomatic-0.39.jar PE " + input1 + " " + input2 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string
            # keep_logging("COMMAND: " + cmdstring, cmdstring, logger, 'debug')
            cmdstring = "trimmomatic PE " + input1 + " " + input2 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string
            keep_logging("COMMAND: " + cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
                print ""
            except sp.CalledProcessError:
                    keep_logging('Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                    sys.exit(1)
            keep_logging('End: Data Pre-processing', 'End: Data Pre-processing', logger, 'info')
        else:
            crop_string = 'CROP:' + crop

            cmdstring = "trimmomatic PE " + input1 + " " + input2 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                    keep_logging('Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                    sys.exit(1)
            keep_logging('End: Data Pre-processing', 'End: Data Pre-processing', logger, 'info')
    else:
        keep_logging('Pre-processing SE reads using Trimmomatic.', 'Pre-processing SE reads using Trimmomatic.', logger, 'info')
        # adapter_file = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("Trimmomatic", Config)['trimmomatic_bin'] + "/" + ConfigSectionMap("Trimmomatic", Config)['adaptor_filepath']
        clean_filenames = out_path + ConfigSectionMap("Trimmomatic", Config)['f_p']
        # changing this parameter for KPC variant analysis for keeping both reads. date: 31 August
        illumina_string = 'ILLUMINACLIP:' + adapter_file + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['seed_mismatches'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['palindrome_clipthreshold'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['simple_clipthreshold']
        sliding_string = 'SLIDINGWINDOW:' + ConfigSectionMap("Trimmomatic", Config)['window_size'] + ConfigSectionMap("Trimmomatic", Config)['colon'] + ConfigSectionMap("Trimmomatic", Config)['window_size_quality']
        minlen_string = 'MINLEN:' + ConfigSectionMap("Trimmomatic", Config)['minlength']
        headcrop_string = 'HEADCROP:' + ConfigSectionMap("Trimmomatic", Config)['headcrop_length']
        if not crop:

            cmdstring = "trimmomatic SE " + input1 + " " + clean_filenames + " " + illumina_string + " " + sliding_string + " " + minlen_string + " " + headcrop_string

            keep_logging("COMMAND: " + cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                keep_logging('Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                sys.exit(1)
            keep_logging('End: Data Pre-processing', 'End: Data Pre-processing', logger, 'info')

        else:
            crop_string = 'CROP:' + crop

            cmdstring = "trimmomatic SE " + input1 + " " + clean_filenames + " " + crop_string + " " + illumina_string + " " + sliding_string + " " + minlen_string

            keep_logging("COMMAND: " + cmdstring, cmdstring, logger, 'debug')
            try:
                call(cmdstring, logger)
            except sp.CalledProcessError:
                    keep_logging('Error in Trimming step. Exiting.', 'Error in Trimming step. Exiting.', logger, 'exception')
                    sys.exit(1)
            keep_logging('End: Data Pre-processing', 'End: Data Pre-processing', logger, 'info')
