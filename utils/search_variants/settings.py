# -*- coding: utf8 -*-

import common

base_quality_offset = 33

min_depth = 10
max_depth = 2500
min_variant_supporting_reads = 4
base_quality_threshold = 15

freq_threshold_for_triallelic_site_checker = 0.03
freq_threshold_for_indel_cover_checker = 0.03

class hetero_germline_variant:
    min_freq = 0.1
    max_freq = 0.9
    min_strand_freq = 0.01
    max_strand_freq = 0.99
    min_avg_base_quality = 20


class somatic_variant:
    min_tumor_freq = 0.05
    sufficient_num_variant_reads = None
    max_normal_freq = 0.1
    max_variant_supporting_reads_normal = 1
    min_strand_freq = 0.1
    max_strand_freq = 0.9
    min_avg_base_quality = 0

debug_number_of_lines = 1000000

def to_str():
    out = u"\tdepth: [{0}, {1}]\n".format(min_depth, max_depth)
    out += u"\tmin variant-supporting reads: {0}\n".format(min_variant_supporting_reads)
    out += u"\tbase quality threshold: {0}\n".format(base_quality_threshold)
    out += u"\tmin freq. for triallelic-site check: {0}\n".format(freq_threshold_for_triallelic_site_checker)
    out += u"\tmin freq. for indel-cover check: {0}\n".format(freq_threshold_for_indel_cover_checker)
    out += u"\n"
    out += u"\tfor hetero germline variants:\n"
    out += u"\t\tallele freq.: [{0}, {1}]\n".format(hetero_germline_variant.min_freq, hetero_germline_variant.max_freq)
    out += u"\t\tstrand freq.: [{0}, {1}]\n".format(hetero_germline_variant.min_strand_freq, hetero_germline_variant.max_strand_freq)
    out += u"\t\tmin average base quality: {0}\n".format(hetero_germline_variant.min_avg_base_quality)
    out += u"\n"
    out += u"\tfor somatic variants:\n"
    out += u"\t\tallele freq. in tumor >=  {0}\n".format(somatic_variant.min_tumor_freq)
    out += u"\t\tsufficient number of variant-supporting reads: {0}\n".format(somatic_variant.sufficient_num_variant_reads)
    out += u"\t\tallele freq. in normal <= {0}\n".format(somatic_variant.max_normal_freq)
    out += u"\t\tmax variant-supporting reads in normal: {0}\n".format(somatic_variant.max_variant_supporting_reads_normal)
    out += u"\t\tstrand freq.: [{0}, {1}]\n".format(somatic_variant.min_strand_freq, somatic_variant.max_strand_freq)
    out += u"\t\tmin average base quality: {0}\n".format(somatic_variant.min_avg_base_quality)
    return out

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--min_depth", type="int", default=10, help="min depth")
parser.add_option("--max_depth", type="int", default=2500, help="max depth")
parser.add_option("--min_tumor_freq", type="float", default=0.05)
parser.add_option("--sufficient_num_variant_reads", type="int", default=None)
parser.add_option("--max_normal_freq", type="float", default=0.1)
parser.add_option("--max_variant_supporting_reads_normal", type="int", default=1)
