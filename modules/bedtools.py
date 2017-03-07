__author__ = 'alipirani'

import os
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *
from modules.log_modules import *

def bedtools(out_sorted_bam, out_path, analysis, logger, Config):
    base_cmd = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("bedtools", Config)['bedtools_bin'] + "/" + ConfigSectionMap("bedtools", Config)['base_cmd']
    cmd = "%s genomecov -ibam %s -d > %s/%s_coverage.bed" % (base_cmd, out_sorted_bam, out_path, analysis)
    keep_logging("COMMAND: " + cmd, cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Bedtools unmapped step. Exiting.', 'Error in Bedtools unmapped step. Exiting.', logger, 'exception')
        sys.exit(1)
    final_coverage_file = "%s/%s_coverage.bed" % (out_path, analysis)
    return final_coverage_file


# def parse_bed_file(final_bed_unmapped_file):
#     unmapped_positions_array = []
#     with open(final_bed_unmapped_file, 'rU') as fp:
#         for line in fp:
#             line_array = line.split('\t')
# 	    lower_index = int(line_array[1]) + 1
#             upper_index = int(line_array[2]) + 1
#             for positions in range(lower_index,upper_index):
#                 unmapped_positions_array.append(positions)
#     only_unmapped_positions_file = final_bed_unmapped_file + "_positions"
#     f1=open(only_unmapped_positions_file, 'w+')
#     for i in unmapped_positions_array:
#         p_string = str(i) + "\n"
#         f1.write(p_string)
#     return only_unmapped_positions_file
