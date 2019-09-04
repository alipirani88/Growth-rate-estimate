__author__ = 'alipirani'
import os
from modules.log_modules import keep_logging
from modules.logging_subprocess import *

###Pending
#################################################### BWA Alignment ############################################
def align_bwa(base_cmd,forward_clean, out_path, reference, split_field, analysis, files_to_delete):
    print "\n################## BWA Alignment ##################\n"
    cmd = "%s mem -M -R %s -t 8 %s %s > %s/%s_aln.sam" % (base_cmd,split_field, reference, forward_clean, out_path, analysis)
    print "\nRunning:\n [%s] \n" % cmd
    os.system(cmd)
    out_sam = "%s/%s_aln.sam" % (out_path, analysis)
    files_to_delete.append(out_sam)
    if not os.path.isfile(out_sam):
        print "Problem in aligning the reads\n"
        exit()
        usage()
    else:
        print "\n################## END: BWA Alignment ##################\n"
        return out_sam
#################################################### END: BWA Alignment #######################################


