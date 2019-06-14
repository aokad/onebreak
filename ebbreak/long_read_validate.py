#! /usr/bin/env python

from __future__ import print_function
import sys, os, tempfile, subprocess, shutil
import pysam

from .pyssw import *

def ssw_check(target, query, score_ratio_thres = 1.4, start_pos_thres = 0.1, end_pos_thres = 0.9):

    nMatch = 2
    nMismatch = 2
    nOpen = 3
    nExt = 1
    sLibPath = ""
 
    lEle = []
    dRc = {} 
    dEle2Int = {}
    dInt2Ele = {}
    # init DNA score matrix
    lEle = ['A', 'C', 'G', 'T', 'N']
    dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'T', 'c':'G', 'g':'C', 't':'A'} 
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in range(nEleNum**2)]
    for i in range(nEleNum-1):
        for j in range(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = nMatch
            else:
                lScore[i*nEleNum+j] = -nMismatch

    # translate score matrix to ctypes
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore

    # set flag
    nFlag = 0
    ssw = ssw_lib.CSsw(sLibPath)
    supporting_reads = []

    # iterate query sequence
    for sQId,sQSeq,sQQual in read(query):

        sQSeq = sQSeq.strip('\n')

        # build query profile
        qNum = to_int(sQSeq, lEle, dEle2Int)
        qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)

        # for reverse complement
        # import pdb; pdb.set_trace()
        sQRcSeq = ''.join([dRc[x] for x in sQSeq[::-1]])
        qRcNum = to_int(sQRcSeq, lEle, dEle2Int)
        qRcProfile = ssw.ssw_init(qRcNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)

        # set mask len
        if len(sQSeq) > 30:
            nMaskLen = int(len(sQSeq) / 2)
        else:
            nMaskLen = 15

        # iter target sequence
        for sRId,sRSeq,_ in read(target):

            rNum = to_int(sRSeq, lEle, dEle2Int)

            # format ofres: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
            res = align_one(ssw, qProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)

            # align rc query
            resRc = align_one(ssw, qRcProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)
    
            # build cigar and trace back path
            strand = 0
            if resRc == None or res[0] > resRc[0]:
                resPrint = res
                strand = 0
                sCigar, sQ, sA, sR = buildPath(sQSeq, sRSeq, res[4], res[2], res[8])
            else:
                resPrint = resRc
                strand = 1
                sCigar, sQ, sA, sR = buildPath(sQRcSeq, sRSeq, resRc[4], resRc[2], resRc[8])

            # import pdb; pdb.set_trace()
            if int(resPrint[0]) > score_ratio_thres * len(sRSeq) and int(resPrint[2]) + 1 < start_pos_thres * len(sRSeq) and int(resPrint[3]) + 1 > end_pos_thres * len(sRSeq):
                supporting_reads.append([sQId, resPrint[0], resPrint[2] + 1, resPrint[3] + 1])

        ssw.init_destroy(qProfile)
        ssw.init_destroy(qRcProfile)
       
    return(supporting_reads)



def long_read_validate_by_alignment(input_file, output_file, bam_file):

    bam_ps = pysam.AlignmentFile(bam_file, "rb")
 
    rname2key = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')               
            if F[0] == "Chr": continue # header
            tchr, tpos, tdir, tjuncseq = F[0], int(F[1]), F[2], F[3]
            key = ','.join([tchr, str(tpos), tdir, tjuncseq])
            for read in bam_ps.fetch(tchr, tpos - 100, tpos + 100):

                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)


    hout = open(output_file + ".tmp3.query_seq.unsorted", 'w')
    for read in bam_ps.fetch():

        flags = format(int(read.flag), "#014b")[:1:-1]
    
        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads             
        if flags[10] == "1": continue

        if read.qname in rname2key:

            for key in rname2key[read.qname]:

                print(key + '\t' + read.qname + '\t' + read.query_sequence, file = hout)

    hout.close()
    bam_ps.close()

    hout = open(output_file + ".tmp3.query_seq.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp3.query_seq.unsorted"], stdout = hout)
    hout.close()

    

    key2contig = {}
    with open(input_file, 'r') as hin:

        header = hin.readline().rstrip('\n').split('\t')
        header2ind = {}
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr, tpos, tdir, tjuncseq = F[0], int(F[1]), F[2], F[3]
            key = ','.join([tchr, str(tpos), tdir, tjuncseq])

            post_contig = F[header2ind["Contig_Post_BP"]]
            all_contig = F[header2ind["Contig_All"]]
            pre_contig = F[header2ind["Contig_All"]][:(len(F[header2ind["Contig_All"]]) - len(F[header2ind["Contig_Post_BP"]]))]

            if len(pre_contig) > len(post_contig):
                proc_contig = pre_contig[-len(post_contig):] + post_contig
            else:
                proc_contig = pre_contig + post_contig[:len(pre_contig)]

            key2contig[key] = proc_contig


    tmp_dir = tempfile.mkdtemp()
    # tmp_dir = "tmp"
    # os.makedirs(tmp_dir)

    temp_key = ""
    temp_id2seq = {}
    temp_junc_seq = ""
    temp_total_read_count = 0
    key2sread_count = {}
    key2sread_count_all = {}
    with open(output_file + ".tmp3.query_seq.sorted") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            # print('>' + F[0] + '\n' + key2contig[F[0]])
            with open(tmp_dir + '/' + F[0] + ".query", 'w') as hout1:
                print('>' + F[0] + '\n' + key2contig[F[0]], file = hout1) 

            if temp_key != F[0]:
                if temp_key != "" and len(key2contig[temp_key]) >= 100:
                    print(temp_key)
                    # print(len(key2contig[temp_key]))
                    supporting_reads = ssw_check(tmp_dir + '/' + temp_key + ".query", tmp_dir + '/' + temp_key + ".target")
                    key2sread_count[temp_key] = len(supporting_reads)
                    key2sread_count_all[temp_key] = temp_total_read_count
                    # tchr, tpos, tdir, tjuncseq = temp_key.split(',')
                    # key2contig[temp_key], key2contig_all[temp_key] = assemble_seq(temp_id2seq, temp_junc_seq, tjuncseq, output_file)

                hout2 = open(tmp_dir + '/' + F[0] + ".target", 'w')
                temp_key = F[0]
                temp_total_read_count = 0
                FF = temp_key.split(',')
            
            print('>' + F[1] + '\n' + F[2], file = hout2)
            temp_total_read_count = temp_total_read_count + 1

        if key2contig[temp_key] != "":
            print(temp_key)
            supporting_reads = ssw_check(tmp_dir + '/' + temp_key + ".query", tmp_dir + '/' + temp_key + ".target")
            # print(supporting_reads)
            key2sread_count[temp_key] = len(supporting_reads)
            key2sread_count_all[temp_key] = temp_total_read_count

    shutil.rmtree(tmp_dir)

    return([key2sread_count, key2sread_count_all])



def add_long_read_validate(input_file, output_file, tumor_bam_file, control_bam_file = None):

    key2sread_count_tumor, key2sread_count_all_tumor = long_read_validate_by_alignment(input_file, output_file + '.tumor', tumor_bam_file)

    if control_bam_file is not None:
        key2sread_count_control, key2sread_count_all_control = long_read_validate_by_alignment(input_file, output_file + '.control', control_bam_file)

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:

        header = hin.readline().rstrip('\n')
        if control_bam_file is not None:
            print(header + '\t' + "Long_Read_Checked_Read_Num_Tumor" + '\t' + "Long_Read_Supporting_Read_Num_Tumor" + '\t' + \
                                  "Long_Read_Checked_Read_Num_Control" + '\t' + "Long_Read_Supporting_Read_Num_Control", file = hout)
        else:
            print(header + '\t' + "Long_Read_Checked_Read_Num_Tumor" + '\t' + "Long_Read_Supporting_Read_Num_Tumor", file = hout)

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])

            sread_count_tumor = key2sread_count_tumor[key] if key in key2sread_count_tumor else 0
            sread_count_all_tumor = key2sread_count_all_tumor[key] if key in key2sread_count_all_tumor else 0

            if control_bam_file is not None:
                sread_count_control = key2sread_count_control[key] if key in key2sread_count_control else 0
                sread_count_all_control = key2sread_count_all_control[key] if key in key2sread_count_all_control else 0

            if control_bam_file is not None:
                print('\t'.join(F) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor) + '\t' + str(sread_count_all_control) + '\t' + str(sread_count_control), file = hout)
            else:
                print('\t'.join(F) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor), file = hout)

    hout.close()
 
    """
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')    
            key = ','.join(F[:4])

            if key not in key2contig: continue
            contig = key2contig[key]
            contig_all = key2contig_all[key]
            if len(contig) < min_contig_length: continue
            # if contig[:8] != F[3][:8]: continue

            
            print('\t'.join(F) + '\t' + contig + '\t' + contig_all, file = hout)

    hout.close()

    # subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.sorted"])
    # subprocess.call(["rm", "-rf", output_file + ".tmp2.contig.unsorted"])


    """

if __name__ == "__main__":

    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    tumor_bam = sys.argv[3]

    generate_contig(input_file, output_file, tumor_bam)

    """
    import sys
    target = sys.argv[1]
    query = sys.argv[2]
    print(main0(target, query))
    """
       
 
