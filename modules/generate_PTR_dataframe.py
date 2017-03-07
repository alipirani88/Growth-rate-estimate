from __future__ import division
__author__ = 'alipirani'
import statistics
import os
import readline
import argparse
from itertools import islice
import subprocess
import numpy as np
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *
from modules.log_modules import *
from modules.generate_R_scripts import *

def generate_moving_sum_results(moving_sum_array, out_path):
    out_file = out_path + "_moving_sum_bins"
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in moving_sum_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')
    peak = max(moving_sum_array)
    through = min([x for x in moving_sum_array if x !=0])
    PTR_moving = peak/through
    out_file_ptr = out_path.replace('.csv', '_PTR.txt')
    with open(out_file_ptr, 'w+') as out:
        out.write(str(out_path)+' moving_sum_array :\t'+str(PTR_moving)+'\n')

def calculate_per_of_reads(median_sliding_window_array, bedfile):
    stats_file = bedfile.replace('_coverage.bed', '_alignment_stats')
    cmd = "grep \'mapped (\' " + stats_file
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out = out.strip()
    out_split = out.split(' ')
    mapped_reads = out_split[0]
    print "The no of mapped reads are: %s" % mapped_reads
    out_file = os.path.dirname(bedfile) + "/" + os.path.basename(bedfile)[0:20] + "_perc_bins.csv"
    with open(out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                perc = (i*100)/int(mapped_reads)
                out.write(str(count)+','+str(perc)+'\n')

def smoothing_1(read_counts, window, out_path, bedfile, logger, Config):
    keep_logging('Generating data frame for Growth Rate(PTR) analysis', 'Generating data frame for Growth Rate(PTR) analysis', logger, 'info')

    ## Step 1
    raw_count_sliding_window_array = []
    zero_bins = 0
    for i in xrange(0, len(read_counts), 100):
        start = i
        end = i + window
        if read_counts[start:end].count(0) == window:
            zero_bins += 1
        raw_count_sliding_window_array.append(read_counts[start:end])
    keep_logging('The length of raw_count_sliding_window_array is {}'.format(len(raw_count_sliding_window_array)), 'The length of raw_count_sliding_window_array is {}'.format(len(raw_count_sliding_window_array)), logger, 'debug')
    keep_logging('The number of bins with no mapped reads: {}'.format(str(zero_bins)), 'The number of bins with no mapped reads: {}'.format(str(zero_bins)), logger, 'debug')

    ## Step 2
    moving_sum_array = []
    for i in raw_count_sliding_window_array:
        moving_sum_array.append(sum(i))
    generate_moving_sum_results(moving_sum_array, out_path)

    ## Step 3
    median_sliding_window_array = []
    for i in xrange(0, len(read_counts), 100):
        start = i
        end = i + window
        if len(moving_sum_array[start:end]) > 5000:
            median_sliding_window_array.append(statistics.median(moving_sum_array[start:end]))

    ## Step 4
    peak = max(median_sliding_window_array)
    through = min([x for x in median_sliding_window_array if x !=0])
    PTR_median = peak/through
    keep_logging('The peak and trough values for median_sliding_window_array is: {};{}'.format(peak, through), 'The peak and trough values for median_sliding_window_array is: {};{}'.format(peak, through), logger, 'info')

    # print "The genomic location for peak values: %s" % (median_sliding_window_array.index(peak))
    # print "The genomic location for trough values: %s" % (median_sliding_window_array.index(through))

    ## Step 5
    out_file_ptr = out_path.replace('.csv', '_PTR.txt')
    with open(out_file_ptr, 'a') as out:
        out.write(str(out_path)+' median_sliding_window_array :\t'+str(PTR_median)+'\n')
        out.write(str(out_path)+' Peak and trough location:\t' + str(median_sliding_window_array.index(peak)) + '\t' + str(median_sliding_window_array.index(through)) + '\n')

    with open(out_path, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in median_sliding_window_array:
                count += 1
                out.write(str(count)+','+str(i)+'\n')

    ## Step 6
    calculate_per_of_reads(median_sliding_window_array, bedfile)

    return PTR_median

def generate_coverage(bedfile, read_counts):
    partitions = [read_counts[i:i+13079] for i in xrange(0, len(read_counts), 13079)]
    partitions = partitions if len(partitions[-1]) == 13079 else partitions[:-1]
    mean_array = []
    for x in partitions:
        mean = statistics.mean(x)
        mean_array.append(mean)
    coverage_out_file = bedfile + "_coverage"
    with open(coverage_out_file, 'w') as out:
        header = "bin,count\n"
        count = 0
        out.write(header)
        for i in mean_array:
            count += 1
            out.write(str(count)+','+str(i)+'\n')

# Main Method
def generate_PTR_dataframe(bedfile, out_path, logger, Config):
    keep_logging('Reading BED file: {} and generating read counts matrix'.format(bedfile), 'Reading BED file: {} and generating read count matrix'.format(bedfile), logger, 'info')
    out_file = out_path + "/" + os.path.basename(bedfile)[0:20] + "_bins.csv"
    window = 10000
    read_counts = []
    with open(bedfile) as fp:
        for line in fp:
            line = line.strip()
            line_split = line.split('\t')
            counts = int(line_split[2])
            read_counts.append(int(counts))
    PTR_median = smoothing_1(read_counts, window, out_file, bedfile, logger, Config)
    perc_bins_matrix = os.path.dirname(bedfile) + "/" + os.path.basename(bedfile)[0:20] + "_perc_bins.csv"
    generate_perc_coverage_graph(perc_bins_matrix, PTR_median)
    #generate_coverage(bedfile, read_counts)