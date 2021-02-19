__author__ = 'alipirani'

import os
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *
from modules.log_modules import *

def bedtools(out_sorted_bam, out_path, analysis, logger, Config):
    #base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bedtools", Config)['bedtools_bin'] + "/" + ConfigSectionMap("bedtools", Config)['base_cmd']
    cmd = "%s genomecov -ibam %s -d > %s/%s_coverage.bed" % (ConfigSectionMap("bedtools", Config)['base_cmd'], out_sorted_bam, out_path, analysis)
    keep_logging("COMMAND: " + cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Bedtools unmapped step. Exiting.', 'Error in Bedtools unmapped step. Exiting.', logger, 'exception')
        sys.exit(1)
    final_coverage_file = "%s/%s_coverage.bed" % (out_path, analysis)
    return final_coverage_file
