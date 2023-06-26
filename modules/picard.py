__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *

def markduplicates(out_sorted_bam, out_path, analysis, files_to_delete, logger, Config):
    # base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("picard", Config)['picard_bin'] + "/" + ConfigSectionMap("picard", Config)['base_cmd']
    base_cmd = ConfigSectionMap("picard", Config)['base_cmd']
    keep_logging('Removing PCR duplicates using PICARD', 'Removing PCR duplicates using PICARD', logger, 'info')

    # cmd = "java -Xmx1G -jar %s MarkDuplicates REMOVE_DUPLICATES=true INPUT=%s OUTPUT=%s/%s_aln_marked.bam METRICS_FILE=%s/%s_markduplicates_metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT" % (base_cmd, out_sorted_bam, out_path, analysis, out_path, analysis)
    cmd = "%s MarkDuplicates REMOVE_DUPLICATES=true INPUT=%s OUTPUT=%s/%s_aln_marked.bam METRICS_FILE=%s/%s_markduplicates_metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT" % (
    base_cmd, out_sorted_bam, out_path, analysis, out_path, analysis)


    keep_logging("COMMAND: " + cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
            keep_logging('Error in Picard Duplicates Removal step. Exiting.', 'Error in Picard Duplicates Removal step. Exiting.', logger, 'exception')
            sys.exit(1)
    out_marked_bam = "%s/%s_aln_marked.bam" % (out_path, analysis)
    #files_to_delete.append(out_marked_bam)
    if not os.path.isfile(out_marked_bam):
        print "Problem in Picard MarkDuplicate Step\n"
        keep_logging('Problem in Picard MarkDuplicate Step', 'Problem in Picard MarkDuplicate Step', logger, 'exception')
        exit()
    else:
        return out_marked_bam

## picard: Prepare Reference Sequence Dictionary for GATK indel realignment
def picard_seqdict(reference_filename, reference):
    dict_name = os.path.splitext(os.path.basename(reference_filename))[0] + ".dict"
    # cmd = "java -jar %s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s/%s" % (base_cmd, reference_filename, ConfigSectionMap(reference, Config)['ref_path'],dict_name)
    cmd = "%s CreateSequenceDictionary REFERENCE=%s OUTPUT=%s/%s" % (
    base_cmd, reference_filename, ConfigSectionMap(reference, Config)['ref_path'], dict_name)
    print "\nRunning:\n [%s] \n" % cmd
    os.system(cmd)

