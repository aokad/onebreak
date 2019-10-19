#! /usr/bin/env python

from __future__ import print_function
import sys, os, math, re, subprocess
import pysam, swalign
import annot_utils.gene, annot_utils.exon

from .filt import get_target_bp
from . import my_seq

def assemble_seq(readid2seq, pre_juncseq, post_juncseq, tmp_file_path, fermi_lite_option = "-l 20", swalign_score_thres = 35):

    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    hout = open(tmp_file_path + ".tmp3.assemble_input.fa", 'w')
    for tid in sorted(readid2seq):
        print('>' + tid, file = hout)
        print(readid2seq[tid], file = hout)
    hout.close()
   
    fermi_lite_options = fermi_lite_option.split(' ')
    hout = open(tmp_file_path + ".tmp3.assemble_output.fq", 'w')
    sret = subprocess.call(["fml-asm"] + fermi_lite_options + [tmp_file_path + ".tmp3.assemble_input.fa"], stdout = hout) 
    hout.close()

    if sret != 0:
        print("fml-asm error, error code: " + str(sret), file = sys.stderr)
        sys.exit()
 
    line_num = 0
    temp_contig = ""
    temp_contig_all = ""
    with open(tmp_file_path + ".tmp3.assemble_output.fq", 'r') as hin:
        for line in hin:
            line_num = line_num + 1
            if line_num % 4 == 2:
                tseq = line.rstrip('\n')

                aln_1 = sw.align(tseq, pre_juncseq)
                if aln_1.score >= swalign_score_thres:
                    ttcontig = tseq[aln_1.r_end:]
                    if ttcontig[:8] == post_juncseq and len(ttcontig) > len(temp_contig): 
                        temp_contig = ttcontig
                        temp_contig_all = tseq

                aln_2 = sw.align(tseq, my_seq.reverse_complement(pre_juncseq))
                if aln_2.score >= swalign_score_thres:
                    ttcontig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                    if ttcontig[:8] == post_juncseq and len(ttcontig) > len(temp_contig): 
                        temp_contig = ttcontig
                        temp_contig_all = my_seq.reverse_complement(tseq)            


    subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.fa"])
    subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_output.fq"])

    if temp_contig == "": temp_contig = "---" 
    if temp_contig_all == "": temp_contig_all = "---"
    return [temp_contig, temp_contig_all]


def generate_contig(input_file, output_file, tumor_bam, reference_genome, min_contig_length, fermi_lite_option, swalign_length = 20, swalign_score = 35):

    tumor_bam_bh = pysam.Samfile(tumor_bam, "rb")
    
    readid2key = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')               
            if F[0] == "Chr": continue

            total_num_tumor, read_ids_tumor, mapping_quals_tumor, clipping_sizes_tumor, alignment_sizes_tumor, juncseq_base_quals_tumor = get_target_bp(F[0], F[1], F[2], F[3], tumor_bam_bh)

            for read_id in read_ids_tumor.split(';'):
                readid2key[re.sub(r'/\d$', '', read_id)] = ','.join(F[:4])

 
    bamfile = pysam.Samfile(tumor_bam, "rb")

    hout = open(output_file + ".tmp2.contig.unsorted", 'w')
    for read in bamfile.fetch():
       
        if read.qname in readid2key:
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            print(readid2key[read.qname] + '\t' + read.qname + ("/1" if flags[6] == "1" else "/2") + '\t' + read.query_sequence, file = hout)

    hout.close()

    hout = open(output_file + ".tmp2.contig.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp2.contig.unsorted"], stdout = hout)
    hout.close()


    temp_key = ""
    temp_id2seq = {}
    temp_junc_seq = ""
    key2contig = {}
    key2contig_all = {}
    with open(output_file + ".tmp2.contig.sorted") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if temp_key != F[0]:
                if len(temp_id2seq) > 0:
                    # print(temp_key)
                    tchr, tpos, tdir, tjuncseq = temp_key.split(',')
                    key2contig[temp_key], key2contig_all[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file, fermi_lite_option, swalign_score)

                temp_key = F[0]
                temp_id2seq = {}
                FF = temp_key.split(',')
                if FF[2] == "+":
                    temp_junc_seq = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - swalign_length, int(FF[1]))
                else:
                    temp_junc_seq = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + swalign_length))

            temp_id2seq[F[1]] = F[2]

        if len(temp_id2seq) > 0:
            tchr, tpos, tdir, tjuncseq = temp_key.split(',') 
            key2contig[temp_key], key2contig_all[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file, fermi_lite_option, swalign_score)


    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            if F[0] == "Chr" and F[1] == "Pos":
                print('\t'.join(F) + '\t' + "Contig_Post_BP" + '\t' + "Contig_All", file = hout)
                continue

            key = ','.join(F[:4])

            # if key not in key2contig: continue
            contig = key2contig[key] if key in key2contig else "---"
            contig_all = key2contig_all[key] if key in key2contig_all else "---"
            # if len(contig) < min_contig_length: continue

            
            print('\t'.join(F) + '\t' + contig + '\t' + contig_all, file = hout)

    hout.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.sorted"])
    subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.unsorted"])


