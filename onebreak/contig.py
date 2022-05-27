#! /usr/bin/env python

from __future__ import print_function
import sys, os, math, re, subprocess
import pysam
import annot_utils.gene, annot_utils.exon

from .filt import get_target_bp
from . import my_seq
import parasail

def assemble_seq(readid2seq, pre_juncseq, post_juncseq, tmp_file_path, fermi_lite_option = "-l 20", parasail_score_thres = 35, temp_key = "", debug = False):

    match = 2
    mismatch = -1
    user_matrix = parasail.matrix_create("ACGT", match, mismatch)

    hout = open(tmp_file_path + ".tmp3.assemble_input.%s.fa" % temp_key, 'w')
    for tid in sorted(readid2seq):
        print('>' + tid, file = hout)
        print(readid2seq[tid], file = hout)
    hout.close()
   
    fermi_lite_options = fermi_lite_option.split(' ')
    hout = open(tmp_file_path + ".tmp3.assemble_output.%s.fq" % temp_key, 'w')
    sret = subprocess.call(["fml-asm"] + fermi_lite_options + [tmp_file_path + ".tmp3.assemble_input.%s.fa" % temp_key], stdout = hout) 
    hout.close()

    if sret != 0:
        print("fml-asm error, error code: " + str(sret), file = sys.stderr)
        sys.exit()

    line_num = 0
    temp_contig = ""
    temp_contig_all = ""
    with open(tmp_file_path + ".tmp3.assemble_output.%s.fq" % temp_key, 'r') as hin:
        for line in hin:
            line_num = line_num + 1
            if line_num % 4 == 2:
                tseq = line.rstrip('\n')

                aln_1 = parasail.ssw(tseq, pre_juncseq, 1, 1, user_matrix)
                if aln_1.score1 >= parasail_score_thres:
                    ttcontig = tseq[aln_1.read_end1+1:]
                    if ttcontig[:8] == post_juncseq and len(ttcontig) > len(temp_contig): 
                        temp_contig = ttcontig
                        temp_contig_all = tseq

                aln_2 = parasail.ssw(tseq, my_seq.reverse_complement(pre_juncseq), 1, 1, user_matrix)
                if aln_2.score1 >= parasail_score_thres:
                    ttcontig = my_seq.reverse_complement(tseq[:aln_2.read_begin1])
                    if ttcontig[:8] == post_juncseq and len(ttcontig) > len(temp_contig): 
                        temp_contig = ttcontig
                        temp_contig_all = my_seq.reverse_complement(tseq)

    if debug == False:
        subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_input.%s.fa" % temp_key])
        subprocess.call(["rm", "-rf", tmp_file_path + ".tmp3.assemble_output.%s.fq" % temp_key])

    if temp_contig == "": temp_contig = "---" 
    if temp_contig_all == "": temp_contig_all = "---"
    return [temp_contig, temp_contig_all]


def generate_contig(input_file, output_file, tumor_bam, reference_genome, min_contig_length, fermi_lite_option, sort_option, swalign_length = 20, swalign_score = 35, debug = False):

    seq_filename, seq_ext = os.path.splitext(tumor_bam)
    if seq_ext == ".cram":
        tumor_bam_bh = pysam.Samfile(tumor_bam, "rc", reference_filename=reference_genome)
    else:
        tumor_bam_bh = pysam.Samfile(tumor_bam, "rb")
    
    readid2key = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "Chr": continue

            total_num_tumor, read_ids_tumor, mapping_quals_tumor, clipping_sizes_tumor, alignment_sizes_tumor, juncseq_base_quals_tumor = get_target_bp(F[0], F[1], F[2], F[3], tumor_bam_bh)

            for read_id in read_ids_tumor.split(';'):
                readid2key[re.sub(r'/\d$', '', read_id)] = ','.join(F[:4])

 
    if seq_ext == ".cram":
        bamfile = pysam.Samfile(tumor_bam, "rc", reference_filename=reference_genome)
    else:
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
    subprocess.call(["sort", "-k1,1"] + sort_option.split(" ") + [output_file + ".tmp2.contig.unsorted"], stdout = hout)
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
                    tchr, tpos, tdir, tjuncseq = temp_key.split(',')
                    key2contig[temp_key], key2contig_all[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file, fermi_lite_option, swalign_score, temp_key, debug)

                temp_key = F[0]
                temp_id2seq = {}
                FF = temp_key.split(',')
                if FF[2] == "+":
                    temp_junc_seq = my_seq.get_seq(reference_genome, FF[0], int(FF[1]) - swalign_length + 1, int(FF[1]))
                else:
                    temp_junc_seq = my_seq.reverse_complement(my_seq.get_seq(reference_genome, FF[0], int(FF[1]), int(FF[1]) + swalign_length + 1))

            temp_id2seq[F[1]] = F[2]

        if len(temp_id2seq) > 0:
            tchr, tpos, tdir, tjuncseq = temp_key.split(',') 
            key2contig[temp_key], key2contig_all[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file, fermi_lite_option, swalign_score, temp_key, debug)


    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            if F[0] == "Chr" and F[1] == "Pos":
                print('\t'.join(F) + '\t' + "Contig_Post_BP" + '\t' + "Contig_All", file = hout)
                continue

            key = ','.join(F[:4])

            contig = key2contig[key] if key in key2contig else "---"
            contig_all = key2contig_all[key] if key in key2contig_all else "---"

            
            print('\t'.join(F) + '\t' + contig + '\t' + contig_all, file = hout)

    hout.close()

    if debug == False:
        subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.sorted"])
        subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.unsorted"])

def filt_contig(input_file, output_file, reference_genome, parasail_refseq_length = 200, parasail_contig_post_bp_length = 25, parasail_score_thres = 45, debug = False):

    match = 2
    mismatch = -1
    user_matrix = parasail.matrix_create("ACGT", match, mismatch)
    
    hout = open(output_file, 'w')
    with open(input_file) as hin:

        header = hin.readline().rstrip('\n').split('\t')
        print('\t'.join(header + ['ScoreIndel']), file = hout)

        header2ind = {}
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            temp_contig_post_bp = F[header2ind['Contig_Post_BP']]

            if temp_contig_post_bp != '---':
                if len(temp_contig_post_bp) > parasail_contig_post_bp_length:
                    temp_contig_post_bp= temp_contig_post_bp[:parasail_contig_post_bp_length]
                temp_ref_seq = my_seq.get_seq(reference_genome, F[0], int(F[1]) - parasail_refseq_length, int(F[1]) + parasail_refseq_length)
                aln_1 = parasail.ssw(temp_ref_seq, temp_contig_post_bp, 1, 1, user_matrix)
                #hout.write(line.rstrip('\n') + "\t" + "\t".join(["%d" % (aln_1.score1), temp_contig_post_bp, temp_ref_seq]) + "\n")
                hout.write(line.rstrip('\n') + "\t" + "\t".join(["%d" % (aln_1.score1)]) + "\n")
            else:
                hout.write(line.rstrip('\n') + "\t" + "\t".join(["0"]) + "\n")
    hout.close()

if __name__ == '__main__':
    pass
    #input_file = "/home/aiokada/sandbox/onebreak/output/ERP001942_0.1.0b9_control_full/ERR188022/ERR188022.onebreak-contig.txt"
    #output_file = "/home/aiokada/sandbox/onebreak/output/ERP001942_0.1.0b10_control_full/ERR188022/ERR188022.onebreak-contig.txt"
    #reference_genome = "/home/aiokada/resources/database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
    #filt_contig(input_file, output_file, reference_genome, parasail_refseq_length = 200, parasail_contig_post_bp_length = 25, parasail_score_thres = 45, debug = False)
