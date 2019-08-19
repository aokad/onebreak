#! /usr/bin/env python

from __future__ import print_function

import sys, os, gzip, subprocess
# from . import parse, filt, contig, classify
from . import parse, filt, contig, classify, long_read_validate

def parse_main(args):

    parse.parse_bp_from_bam(args.bam_file, args.output_file + ".bp.tmp.txt", args.key_seq_size, args.min_major_clip_size, args.max_minor_clip_size)

    parse.cluster_breakpoint(args.output_file + ".bp.tmp.txt", args.output_file + ".bp.clustered.tmp.txt", args.check_interval)

    hout = open(args.output_file + ".bp.clustered.sorted.tmp.txt", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", args.output_file + ".bp.clustered.tmp.txt"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print("Error in sorting merged junction file", file = sys.stderr)
        sys.exit(1)

    hout = open(args.output_file, 'w')
    s_ret = subprocess.call(["bgzip", "-f", "-c", args.output_file + ".bp.clustered.sorted.tmp.txt"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print("Error in compression merged junction file", file = sys.stderr)
        sys.exit(1)


    s_ret = subprocess.call(["tabix", "-p", "bed", args.output_file])
    if s_ret != 0:
        print("Error in indexing merged junction file", file = sys.stderr)
        sys.exit(1)

    subprocess.call(["rm", "-f", args.output_file + ".bp.tmp.txt"])
    subprocess.call(["rm", "-f", args.output_file + ".bp.clustered.tmp.txt"])
    subprocess.call(["rm", "-f", args.output_file + ".bp.clustered.sorted.tmp.txt"])
   

def merge_control_main(args):
    
    """
    function for merging control samples
    """
    # make directory for output if necessary
    if os.path.dirname(args.output_file) != "" and not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    # merge control samples
    hout = open(args.output_file + ".unsorted.txt", 'w')
    with open(args.bp_file_list, 'r') as hin:
        for line in hin:
            bp_file = line.rstrip('\n')
            with gzip.open(bp_file, 'r') as hin2:
                for line2 in hin2:
                    F = line2.rstrip('\n').split('\t')
                    support_num = len(F[7].split(';'))
                    if support_num < args.support_num_thres: continue 
                    print(F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3] + '\t' + F[4] + '\t' + str(support_num), file = hout)
                

    hout = open(args.output_file + ".sorted.txt", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", args.output_file + ".unsorted.txt"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print("Error in sorting merged junction file", file = sys.stderr)
        sys.exit(1)


    hout = open(args.output_file + ".merged.txt", 'w')
    with open(args.output_file + ".sorted.txt", 'r') as hin:
        temp_key = ""
        temp_read_num = []
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3]
            support_num = int(F[5])
            if key != temp_key:
                if temp_key != "":
                    if len(temp_read_num) >= args.sample_num_thres:
                        print(temp_key + '\t' + ','.join(temp_read_num), file = hout)
                temp_key = key
                temp_read_num = []
            temp_read_num.append(str(support_num))

        if len(temp_read_num) >= args.sample_num_thres:
            print(temp_key + '\t' + ','.join(temp_read_num), file = hout)

    hout.close()



    hout = open(args.output_file, 'w')
    s_ret = subprocess.call(["bgzip", "-f", "-c", args.output_file + ".merged.txt"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print("Error in compression merged junction file", file = sys.stderr)
        sys.exit(1)


    s_ret = subprocess.call(["tabix", "-p", "bed", args.output_file])
    if s_ret != 0:
        print("Error in indexing merged junction file", file = sys.stderr)
        sys.exit(1)

    subprocess.call(["rm", "-f", args.output_file + ".unsorted.txt"])
    subprocess.call(["rm", "-f", args.output_file + ".sorted.txt"])
    subprocess.call(["rm", "-f", args.output_file + ".merged.txt"])


def filt_main(args):

    filt.filter_by_merged_control(args.tumor_bp_file, args.output_file + ".tmp.filt1.txt", args.merged_control_file,
                                  args.min_median_mapq, args.min_max_clip_size, args.min_second_juncseq_baseq, args.permissible_range)

    filt.filter_by_peakedness(args.output_file + ".tmp.filt1.txt", args.tumor_bp_file, args.output_file + ".tmp.filt2.txt")

    filt.filter_by_matched_control_local_check(args.output_file + ".tmp.filt2.txt", args.matched_control_bp_file, args.output_file + ".tmp.filt3.txt")

    filt.filter_by_allele_freq(args.output_file + ".tmp.filt3.txt", args.output_file, 
                               args.tumor_bam, args.matched_control_bam, args.min_variant_num_tumor, args.min_VAF_tumor, 
                               args.max_variant_num_control, args.max_VAF_control, args.max_fisher_pvalue)

    if not args.debug:
        subprocess.call(["rm", args.output_file + ".tmp.filt1.txt"])
        subprocess.call(["rm", args.output_file + ".tmp.filt2.txt"])
        subprocess.call(["rm", args.output_file + ".tmp.filt3.txt"])


def long_read_validate_main(args):

    long_read_validate.add_long_read_validate(args.contig_result_file, args.output_file, args.reference_genome, args.tumor_bam, args.control_bam)


def contig_main(args):

    contig.generate_contig(args.tumor_bp_filt_file, args.output_file,
                           args.tumor_bam, args.reference_genome, args.min_contig_length)

    # contig.alignment_contig(args.tumor_bp_filt_file, args.output_file + ".tmp.filt4.txt", args.output_file, 
    #                         args.reference_genome, args.blat_option, args.blat_ooc, args.virus_db, args.repeat_db, args.mitochondria_db, args.adapter_db)
    


    # subprocess.call(["rm", "-f", args.output_file + ".tmp.filt4.txt"])



def classify_main(args):

    classify.classify_by_contig_alignment(args.contig_result_file, args.output_file, args.reference_genome, args.repeat_seq)

    # classify.classify_canonicalSV(args.tumor_bp_contig_file, args.output_canonical_file + ".tmp1.txt")
    
    # classify.filter_doublecount(args.output_canonical_file + ".tmp1.txt", args.output_canonical_file + ".tmp2.txt")

    # classify.filter_rna_junction(args.output_canonical_file + ".tmp2.txt", args.output_canonical_file + ".tmp3.txt", args.grc, args.genome_id)

    # classify.annot_canonicalSV(args.output_canonical_file + ".tmp3.txt", args.output_canonical_file, args.grc, args.genome_id)

    # classify.classify_non_canonicalSV(args.tumor_bp_contig_file, args.output_non_canonical_file + ".tmp1.txt", args.output_canonical_file + ".tmp1.txt")

    # classify.annotate_break_point(args.output_non_canonical_file + ".tmp1.txt", args.output_non_canonical_file + ".tmp2.txt", args.grc, args.genome_id)

    #subprocess.call(["rm", "-f", args.output_canonical_file + ".tmp1.txt"])
    #subprocess.call(["rm", "-f", args.output_canonical_file + ".tmp2.txt"])
