#! /usr/bin/env python

from __future__ import print_function
import sys, os, math, re, subprocess
import pysam, swalign
import annot_utils.gene, annot_utils.exon

from .filt import get_target_bp
from . import my_seq

def assemble_seq(readid2seq, pre_juncseq, post_juncseq, tmp_file_path):

    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)

    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

    hout = open(tmp_file_path + ".tmp3.assemble_input.fa", 'w')
    for tid in sorted(readid2seq):
        print('>' + tid, file = hout)
        print(readid2seq[tid], file = hout)
    hout.close()
    
    hout = open(tmp_file_path + ".tmp3.assemble_output.fq", 'w')
    sret = subprocess.call(["fml-asm", '-l', '20', tmp_file_path + ".tmp3.assemble_input.fa"], stdout = hout) 
    hout.close()

    if sret != 0:
        print("fml-asm error, error code: " + str(sret), file = sys.stderr)
        sys.exit()
 
    line_num = 0
    temp_contig = ""
    with open(tmp_file_path + ".tmp3.assemble_output.fq", 'r') as hin:
        for line in hin:
            line_num = line_num + 1
            if line_num % 4 == 2:
                tseq = line.rstrip('\n')

                aln_1 = sw.align(tseq, pre_juncseq)
                if aln_1.score >= 35:
                    ttcontig = tseq[aln_1.r_end:]
                    if ttcontig[:8] == post_juncseq and len(ttcontig) > len(temp_contig): temp_contig = ttcontig
                
                aln_2 = sw.align(tseq, my_seq.reverse_complement(pre_juncseq))
                if aln_2.score >= 35:
                    ttcontig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                    if ttcontig[:8] == post_juncseq and len(ttcontig) > len(temp_contig): temp_contig = ttcontig

    """
    line_num = 0
    temp_score = 35
    temp_contig = ""
    # temp_all_contig = ""
    # temp_genome_contig = ""
    with open(tmp_file_path + ".tmp3.assemble_output.fq", 'r') as hin:
        for line in hin:
            line_num = line_num + 1
            if line_num % 4 == 2:
                tseq = line.rstrip('\n')
                
                aln_1 = sw.align(tseq, junc_seq)
                aln_2 = sw.align(tseq, my_seq.reverse_complement(junc_seq))
            
                if aln_1.score >= aln_2.score:
                    if aln_1.score > temp_score:
                        temp_contig = tseq[aln_1.r_end:]
                        # temp_all_contig = str(tseq)
                        # temp_genome_contig = tseq[:aln_1.r_end]
                        temp_score = aln_1.score

                    elif aln_1.score == temp_score:
                        if len(tseq[aln_1.r_end:]) > len(temp_contig): 
                            temp_contig = tseq[aln_1.r_end:]
                            # temp_all_contig = str(tseq)
                            # temp_genome_contig = tseq[:aln_1.r_end]
                            temp_score = aln_1.score

                elif aln_2.score > aln_1.score:
                    if aln_2.score > temp_score:
                        temp_contig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                        # temp_all_contig = my_seq.reverse_complement(tseq)
                        # temp_genome_contig = my_seq.reverse_complement(tseq[aln_2.r_pos:])
                        temp_score = aln_2.score

                    if aln_2.score == temp_score:
                        if len(tseq[aln_2.r_end:]) > len(temp_contig): 
                            temp_contig = my_seq.reverse_complement(tseq[:aln_2.r_pos])
                            # temp_all_contig = my_seq.reverse_complement(tseq)
                            # temp_genome_contig = my_seq.reverse_complement(tseq[aln_2.r_pos:])
                            temp_score = aln_2.score
    """

    subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.fa"])
    subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_output.fq"])

    return temp_contig


def generate_contig(input_file, output_file, tumor_bp_file, tumor_bam, reference_genome, min_contig_length, swalign_length, swalign_score):
# def generate_contig(input_file, output_file, tumor_bam, reference_genome, min_contig_length):

    """
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
    """

    temp_key = ""
    temp_id2seq = {}
    temp_junc_seq = ""
    key2contig = {}
    with open(output_file + ".tmp2.contig.sorted") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if temp_key != F[0]:
                if len(temp_id2seq) > 0:
                    print(temp_key)
                    tchr, tpos, tdir, tjuncseq = temp_key.split(',')
                    key2contig[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file)

                temp_key = F[0]
                temp_id2seq = {}
                FF = temp_key.split(',')
                if FF[2] == "+":
                    temp_junc_seq = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - 20, int(FF[1]))
                else:
                    temp_junc_seq = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + 20))

            temp_id2seq[F[1]] = F[2]

        if len(temp_id2seq) > 0:
            tchr, tpos, tdir, tjuncseq = temp_key.split(',') 
            key2contig[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file)


    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            key = ','.join(F[:4])

            if key not in key2contig: continue
            contig = key2contig[key]
            if len(contig) < min_contig_length: continue
            # if contig[:8] != F[3][:8]: continue

            
            print('\t'.join(F) + '\t' + contig, file = hout)

    hout.close()

    # subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.sorted"])
    # subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.unsorted"])



