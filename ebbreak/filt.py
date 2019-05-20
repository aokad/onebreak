#! /usr/bin/env python

from __future__ import print_function

import sys, gzip, math, numpy
import pysam
from scipy import stats
from statistics import mean
from . import my_seq


def get_target_bp(target_chr, target_pos, target_dir, target_junc_seq, bamfile):

    total_read_num = 0
    read_ids = []
    mapping_quals = []  
    clipping_sizes = []
    alignment_sizes = []    
    juncseq_base_quals = []

    key_seq_size = len(target_junc_seq)

    # maybe add the regional extraction of bam files
    target_region = target_chr + ':' + target_pos + '-' + target_pos
    for read in bamfile.fetch(target_chr, int(target_pos) - 1, int(target_pos)):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if not aligned
        if flags[2] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # no clipping
        if len(read.cigar) == 1: continue

        total_read_num = total_read_num + 1

        # get the clipping size in the both side
        left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
        right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

        if target_dir == '+' and right_clipping < 2: continue
        if target_dir == '-' and left_clipping < 2: continue

        # if left_clipping < min_major_clip_size and right_clipping < min_major_clip_size: continue

        # get the alignment basic information
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")

        # when the right side is clipped...
        if target_dir == '+':
            clipLen_current = right_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
            juncChr_current = chr_current
            juncPos_current = pos_current + alignmentSize_current - 1

            juncDir_current = "+"
            juncseq_start = readLength_current - clipLen_current
            juncseq_end = readLength_current - max(clipLen_current - key_seq_size, 0)
            juncseq = read.seq[juncseq_start:juncseq_end]
            juncseq_baseq = mean(read.query_qualities[juncseq_start:juncseq_end])

            if str(juncPos_current) != str(target_pos): continue
            if juncseq != target_junc_seq[:len(juncseq)]: continue

            read_ids.append(read.qname + ("/1" if flags[6] == "1" else "/2"))
            mapping_quals.append(str(read.mapq))
            clipping_sizes.append(str(right_clipping))
            alignment_sizes.append(str(alignmentSize_current))
            juncseq_base_quals.append(str(juncseq_baseq))


        if target_dir == '-':

            clipLen_current = left_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
     
            juncChr_current = chr_current
            juncPos_current = pos_current
            juncDir_current = "-"

            juncseq_end = clipLen_current
            juncseq_start = max(clipLen_current - key_seq_size, 0)
            juncseq = my_seq.reverse_complement(read.seq[juncseq_start:juncseq_end])
            juncseq_baseq = mean(read.query_qualities[juncseq_start:juncseq_end])

            if str(juncPos_current) != str(target_pos): continue
            if juncseq != target_junc_seq[:len(juncseq)]: continue

            read_ids.append(read.qname + ("/1" if flags[6] == "1" else "/2"))
            mapping_quals.append(str(read.mapq))
            clipping_sizes.append(str(left_clipping))
            alignment_sizes.append(str(alignmentSize_current))
            juncseq_base_quals.append(str(juncseq_baseq))


    return([total_read_num, ';'.join(read_ids), ';'.join(mapping_quals), ';'.join(clipping_sizes), ';'.join(alignment_sizes), ';'.join(juncseq_base_quals)])




def filter_by_merged_control(tumor_bp_file, output_file, merged_control_file,
                             min_median_mapq, min_max_clip_size, min_second_juncseq_baseq, permissible_range):

    """
    filter by base quality (1st step)
    filtering with merged control
    """

    use_merged_control = True if merged_control_file != "" else False
    if use_merged_control: merged_control_db = pysam.TabixFile(merged_control_file)

    hout = open(output_file, 'w')
    with gzip.open(tumor_bp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mapqs = [int(x) for x in F[6].split(';')]
            clip_sizes = [int(x) for x in F[7].split(';')]
            base_qualities = [float(x) for x in F[9].split(';')]

            # remove breakpoint if supporting read does not meet the criteria below
            # if len(mapqs) == 1: continue
            if len(list(set(clip_sizes))) == 1: continue
            if numpy.median(mapqs) < min_median_mapq: continue
            if max(clip_sizes) < min_max_clip_size: continue
            #if numpy.median(base_qualities) < 20: continue
            if numpy.sort(base_qualities)[-2] < min_second_juncseq_baseq: continue

            # filtering using merged control file
            merged_control_filt_flag = False 
            if use_merged_control:
                tabixErrorFlag = 0
                try:
                    # records = merged_control_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
                    start_pos = int(F[1]) - 1 - permissible_range
                    if start_pos < 0: start_pos = 0
                    records = merged_control_db.fetch(F[0], start_pos , int(F[1]) + 1 + permissible_range)
                except Exception as inst:
                    print("%s: %s" % (type(inst), inst.args), file = sys.stderr)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        # if record[0] == F[0] and record[1] == F[1] and record[2] == F[2]:
                        if record[0] == F[0] and (int(F[1]) - permissible_range) <= int(record[1]) and int(record[1]) <= (int(F[1]) + permissible_range) and record[3] == F[3]:
                            #if ignore_seq_consist or record[4] == F[4]:
                            merged_control_filt_flag = True

            if merged_control_filt_flag: continue

            # print(F[0] + '\t' + str(int(F[1])+1) + '\t' + F[3] + '\t' + F[4], file = hout)
            print('\t'.join(F), file = hout)

    hout.close()



def filter_by_allele_freq(input_file, output_file, tumor_bam, matched_control_bam, min_variant_num_tumor, min_VAF_tumor,
                          max_variant_num_control, max_VAF_control, max_fisher_pvalue):

    """
    filtering by allele frequency
    """

    hout = open(output_file, 'w')

    print('\t'.join(["Chr", "Pos", "Dir", "Junc_Seq", 
                     "Num_Tumor_Total_Read", "Num_Tumor_Var_Read", "Num_Control_Total_Read", "Num_Control_Var_Read",
                     "Minus_Log_Fisher_P_value"]), file = hout)

    tumor_bam_bh = pysam.Samfile(tumor_bam, "rb")
    matched_control_bam_bh = pysam.Samfile(matched_control_bam, "rb")


    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            total_num_tumor, read_ids_tumor, mapping_quals_tumor, clipping_sizes_tumor, alignment_sizes_tumor, juncseq_base_quals_tumor = get_target_bp(F[0], F[2], F[3], F[4], tumor_bam_bh)
            variant_num_tumor = len(read_ids_tumor.split(';'))
            if variant_num_tumor < min_variant_num_tumor: continue    

            if total_num_tumor == 0:
                print("total_num_tumor is 0", file = sys.stderr)
                sys.exit(1)
        
    
            VAF_tumor = float(variant_num_tumor) / total_num_tumor
            if VAF_tumor < min_VAF_tumor: continue


            if matched_control_bam != "":
                total_num_control, read_ids_control, mapping_quals_control, clipping_sizes_control, alignment_sizes_control, juncseq_base_quals_control = get_target_bp(F[0], F[2], F[3], F[4], matched_control_bam_bh)
                variant_num_control = len(read_ids_control.split(';')) if read_ids_control != '' else 0
                VAF_control = float(variant_num_control) / total_num_control if total_num_control > 0 else 1.0

            else:
                variant_read_num_control = "---"
                total_read_num_control = "---"
                VAF_control = "---"

            if variant_num_control > max_variant_num_control: continue
            if VAF_control != "---" and VAF_control > max_VAF_control: continue 

            lpvalue = "---"
            if VAF_control != "":
                oddsratio, pvalue = stats.fisher_exact([[total_num_tumor - variant_num_tumor, variant_num_tumor], [total_num_control - variant_num_control, variant_num_control]], 'less')
                if pvalue < 1e-100: pvalue = 1e-100
                lpvalue = (- math.log(pvalue, 10) if pvalue < 1 else 0)
                lpvalue = round(lpvalue, 4) 

                if 10**(-lpvalue) > float(max_fisher_pvalue): continue

            if VAF_tumor != "---": VAF_tumor = str(round(VAF_tumor, 4))
            if VAF_control != "---": VAF_control = str(round(VAF_control, 4))

            print('\t'.join(F[:5]) + '\t' + str(total_num_tumor) + '\t' + str(variant_num_tumor) + '\t' + VAF_tumor + '\t' + \
                  str(total_num_control) + '\t' + str(variant_num_control) + '\t' + VAF_control + '\t' + str(lpvalue), file = hout)

    hout.close()



def filter_by_base_quality(input_file, output_file, tumor_bam, min_support_num, permissible_range):

    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            tumor_bamfile = pysam.Samfile(tumor_bam, "rb")
            key_baseq = []
            tumor_var_read = 0
            tabixErrorFlag = 0
            try:
                records = tumor_bamfile.fetch(F[0], max(int(F[1]) - 1, 0) , int(F[1]) + 1)
            
            except Exception as inst:
                print("%s: %s" % (type(inst), inst.args), file = sys.stderr)
                tabixErrorMsg = str(inst.args)
                tabixErrorFlag = 1

            if tabixErrorFlag == 0:
                for read in records:
                    flags = format(int(read.flag), "#014b")[:1:-1]

                    # skip if not aligned
                    if flags[2] == "1": continue

                    # skip supplementary alignment
                    if flags[8] == "1" or flags[11] == "1": continue

                    # skip duplicated reads
                    if flags[10] == "1": continue

                    # no clipping
                    if len(read.cigar) == 1: continue

                    left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
                    right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)
                    #print right_clipping

                    if F[2] == "+":
                        if right_clipping < 2: continue

                        juncPos_current = int(read.pos + 1) + read.alen - 1

                        if int(juncPos_current) != int(F[1]): continue

                        juncseq_start = read.rlen - right_clipping
                        juncseq_end = min(juncseq_start + 8, read.rlen)
                        juncseq = read.seq[juncseq_start:juncseq_end]

                        #if numpy.mean(read.query_qualities[juncseq_start:juncseq_end]) < 10:
                        #    continue

                        if juncseq in F[3]: 
                            tumor_var_read += 1
                            key_baseq.append(str(numpy.mean(read.query_qualities[juncseq_start:juncseq_end])))

                    if F[2] == "-":
                        if left_clipping < 2: continue
                        
                        juncPos_current = int(read.pos + 1)

                        if int(juncPos_current) != int(F[1]): continue

                        juncseq_start = left_clipping
                        juncseq_end = max(left_clipping - 8, 0)
                        juncseq = my_seq.reverse_complement(read.seq[juncseq_start:juncseq_end])
                        
                        #if numpy.mean(read.query_qualities[juncseq_start:juncseq_end]) < 10:
                        #    continue
                        
                        if juncseq in F[3]: 
                            tumor_var_read += 1
                            key_baseq.append(str(numpy.mean(read.query_qualities[juncseq_end:juncseq_start])))


            if tumor_var_read < min_support_num: continue
            if numpy.median(list(map(float, key_baseq))) < 20: continue
            if numpy.sort(list(map(float, key_baseq)))[-2] < 30: continue

            print(F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3]  + '\t' + str(tumor_var_read), file = hout)

    hin.close()
    hout.close()



def filter_by_matched_control(input_file, output_file, matched_control_bp_file, normal_bam, max_control_num_thres, permissible_range):

    """
    filtering with matched control
    """

    use_matched_control = True if matched_control_bp_file != "" else False
    if use_matched_control: normal_bamfile = pysam.Samfile(normal_bam, "rb")


    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        for line in hin:
            
            F = line.rstrip('\n').split('\t')
            normal_var_read = 0

            if use_matched_control:

                tabixErrorFlag = 0
                try:
                    records = normal_bamfile.fetch(F[0], max(int(F[1]) - 1, 0) , int(F[1]) + 1)
                
                except Exception as inst:
                    print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

                if tabixErrorFlag == 0:
                    for read in records:
                            flags = format(int(read.flag), "#014b")[:1:-1]

                            # skip if not aligned
                            if flags[2] == "1": continue

                            # skip supplementary alignment
                            if flags[8] == "1" or flags[11] == "1": continue

                            # skip duplicated reads
                            if flags[10] == "1": continue

                            # no clipping
                            if len(read.cigar) == 1: continue

                            if F[2] == "+":
                                right_clipping2 = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)
                                
                                if right_clipping2 < 2: continue

                                juncPos_current2 = int(read.pos + 1) + read.alen - 1
                                if max(int(F[1]) - permissible_range, 0) <= int(juncPos_current2) and int(juncPos_current2) <= (int(F[1]) + permissible_range):

                                    #if juncseq2 in F[3]: 
                                    normal_var_read += 1

                            if F[2] == "-":
                                left_clipping2 = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)

                                if left_clipping2 < 2: continue
                                juncPos_current2 = int(read.pos + 1)
                                if max(int(F[1]) - permissible_range, 0) <= int(juncPos_current2) and int(juncPos_current2) <= (int(F[1]) + permissible_range):

                                   #if juncseq2 in F[3]: 
                                    normal_var_read += 1

            else:
                normal_var_read == "---"

            if use_matched_control and normal_var_read > max_control_num_thres: continue

            print(F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[3]  + '\t' + F[4] + '\t' + str(normal_var_read), file = hout)


"""
def get_target_bp(target_chr, target_pos, target_dir, target_junc_seq, bamfile):

    total_read_num = 0
    read_ids = []
    mapping_quals = []  
    clipping_sizes = []
    alignment_sizes = []    
    juncseq_base_quals = []

    key_seq_size = len(target_junc_seq)

    # maybe add the regional extraction of bam files
    target_region = target_chr + ':' + target_pos + '-' + target_pos
    for read in bamfile.fetch(target_chr, int(target_pos) - 1, int(target_pos)):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if not aligned
        if flags[2] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        # no clipping
        if len(read.cigar) == 1: continue

        total_read_num = total_read_num + 1

        # get the clipping size in the both side
        left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
        right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

        if target_dir == '+' and right_clipping < 2: continue
        if target_dir == '-' and left_clipping < 2: continue

        # if left_clipping < min_major_clip_size and right_clipping < min_major_clip_size: continue

        # get the alignment basic information
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")

        # when the right side is clipped...
        if target_dir == '+':
        # if right_clipping >= min_major_clip_size:
            clipLen_current = right_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
            juncChr_current = chr_current
            juncPos_current = pos_current + alignmentSize_current - 1

            juncDir_current = "+"
            juncseq_start = readLength_current - clipLen_current
            juncseq_end = readLength_current - max(clipLen_current - key_seq_size, 0)
            juncseq = read.seq[juncseq_start:juncseq_end]
            juncseq_baseq = mean(read.query_qualities[juncseq_start:juncseq_end])

            if str(juncPos_current) != str(target_pos): continue

            if len(juncseq) < 8:
                import pdb; pdb.set_trace()

            if juncseq != target_junc_seq[:len(juncseq)]: continue

            read_ids.append(read.qname + ("/1" if flags[6] == "1" else "/2"))
            mapping_quals.append(str(read.mapq))
            clipping_sizes.append(str(right_clipping))
            alignment_sizes.append(str(alignmentSize_current))
            juncseq_base_quals.append(str(juncseq_baseq))

            # read_infos.append('\t'.join([read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), str(right_clipping), str(alignmentSize_current), str(juncseq_baseq)]))
            # print('\t'.join([juncChr_current, str(juncPos_current-1), str(juncPos_current), juncDir_current, juncseq, 
            #       read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), str(right_clipping), str(alignmentSize_current), str(juncseq_baseq)]), file = hout)

        if target_dir == '-':
        # if left_clipping >= min_major_clip_size:

            clipLen_current = left_clipping
            alignmentSize_current = read.alen
            readLength_current = read.rlen
     
            juncChr_current = chr_current
            juncPos_current = pos_current
            juncDir_current = "-"

            juncseq_end = clipLen_current
            juncseq_start = max(clipLen_current - key_seq_size, 0)
            juncseq = my_seq.reverse_complement(read.seq[juncseq_start:juncseq_end])
            juncseq_baseq = mean(read.query_qualities[juncseq_start:juncseq_end])

            if str(juncPos_current) != str(target_pos): continue
            if len(juncseq) < 8:
                import pdb; pdb.set_trace()
            if juncseq != target_junc_seq[:len(juncseq)]: continue

            read_ids.append(read.qname + ("/1" if flags[6] == "1" else "/2"))
            mapping_quals.append(str(read.mapq))
            clipping_sizes.append(str(left_clipping))
            alignment_sizes.append(str(alignmentSize_current))
            juncseq_base_quals.append(str(juncseq_baseq))

            # read_infos.append('\t'.join([read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), str(left_clipping), str(alignmentSize_current), str(juncseq_baseq)]))
            # print('\t'.join([juncChr_current, str(juncPos_current-1), str(juncPos_current), juncDir_current, juncseq, 
            #       read.qname + ("/1" if flags[6] == "1" else "/2"), str(read.mapq), str(left_clipping), str(alignmentSize_current), str(juncseq_baseq)]), file = hout)


    return([total_read_num, ';'.join(read_ids), ';'.join(mapping_quals), ';'.join(clipping_sizes), ';'.join(alignment_sizes), ';'.join(juncseq_base_quals)])
"""
