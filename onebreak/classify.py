#! /usr/bin/env python

import sys, os, re, copy, subprocess, tempfile, shutil
import pysam

from .check_ref import check_reference
import annot_utils


def get_margin(read):

    def get_match_score(cigar_string, nm):
        match_score = 0
        cigar_matches = re.findall(r'\d+[MIDSHN]', cigar_string)
        for cigar_match in cigar_matches:
            if cigar_match[-1] == 'M':
                match_score = match_score + int(cigar_match[:-1])
            elif cigar_match[-1] == 'D':
                match_score = match_score - int(cigar_match[:-1])
            elif cigar_match[-1] == 'I':
                match_score = match_score - int(cigar_match[:-1])

        match_score = match_score - nm

        return(match_score)

    # return 0 for reads with 0 mapping qualities
    if read.mapping_quality == 0:
        return(0)

    margin = float("Inf")
    if read.has_tag("XA"):
    
        target_nm = read.get_tag("NM") if read.has_tag("NM") else 0
        target_match_score = get_match_score(read.cigarstring, target_nm)
 
        XA_tags = read.get_tag("XA").split(';')
        for XA_tag in XA_tags:
            if XA_tag == '': continue
            xchr, xpos, xcigar, xnm = XA_tag.split(',')        
            match_score = get_match_score(xcigar, int(xnm))
            tmargin = max(target_match_score - match_score, 0)
            if tmargin < margin:
                margin = tmargin        

    return(margin)


def get_query_pos(read):

    query_len = 0
    for elm in read.cigartuples:
        if elm[0] in [0, 1, 4, 5]: query_len = query_len + elm[1]

    if read.cigartuples[0][0] == 0:
        query_start = 1
    elif read.cigartuples[0][0] in [4, 5] and read.cigartuples[1][0] == 0:
        query_start = read.cigartuples[0][1] + 1
    else:
        print("cigar is strange", file = sys.stderr)
        sys.exit(1)

    if read.cigartuples[-1][0] == 0:
        query_end = query_len
    elif read.cigartuples[-1][0] in [4, 5] and read.cigartuples[-2][0] == 0:
        query_end = query_len - read.cigartuples[-1][1]
    else:
        print("cigar is strange", file = sys.stderr)
        sys.exit(1)

    if not read.is_reverse:
        return(query_start, query_end)
    else:
        return(query_len - query_end + 1, query_len - query_start + 1)

 
def get_alignment_pos(read_set):

    # no aligned read
    aligned_flag = False
    for read in read_set:
        if not read.is_unmapped: aligned_flag = True

    if aligned_flag == False:
        return([])

    bp_pos2alignment = [] 
    # first check primary alignment
    for read in read_set:
        # if not read.is_secondary and not read.is_supplementary:
        is_primary = 'p' if not read.is_secondary and not read.is_supplementary else 's'
        mapq = read.mapping_quality

        bp_start, bp_end = get_query_pos(read)

        alignment = read.reference_name + ',' + ('+' if read.is_reverse == False else '-') + \
                    str(read.reference_start + 1) + '-' + str(read.reference_end) + \
                    ',' + str(mapq) + ',' + is_primary

        bp_pos2alignment.append([bp_start, bp_end, alignment])


    return(sorted(bp_pos2alignment, key = lambda x: x[1]))



def get_bp_type(read_set):

    # no aligned read
    aligned_flag = False
    for read in read_set:
        if not read.is_unmapped: aligned_flag = True
    
    if aligned_flag == False:
        return("Unmapped")

    bp_type = "Unclassified"

    # first check primary alignment margin
    for read in read_set:
        if not read.is_secondary and not read.is_supplementary:
            primary_margin = get_margin(read)


    if primary_margin > 5:

        bp_type = "Unique_Mapped"
        if len(read_set) == 1:

            # plain SV
            left_clipping_size = 0
            match = re.search(r'^(\d+)S', read.cigarstring)
            if match is not None:
                left_clipping_size = int(match.group(1))
        
            right_clipping_size = 0
            match = re.search(r'(\d+)S$', read.cigarstring)
            if match is not None:
                right_clipping_size = int(match.group(1))

            if not read.is_reverse and left_clipping_size == 0 and right_clipping_size < 5: bp_type = "Plain SV"
            if read.is_reverse and left_clipping_size < 5 and right_clipping_size == 0: bp_type = "Plain SV"

            if not read.is_reverse and left_clipping_size > 0 and right_clipping_size < 5: bp_type = "Inseq SV"
            if read.is_reverse and left_clipping_size < 5 and right_clipping_size > 0: bp_type = "Inseq SV"

        else:

            bp_type = "Complex"

    else:
        bp_type = "Noncanonical"

    return(bp_type)


def get_sv_info(bp_key, bp_start, first_alignment):

    bchr, bpos, bdir, bseq = bp_key.split(',')

    tchr, tregion, _, _ = first_alignment.split(',')
    tdir = tregion[0]
    tstart, tend = tregion[1:].split('-')
 
    op_pos = tstart if tdir == '+' else tend

    
    if bchr != tchr:
        sv_type = "Translocation"
        sv_size = '---'
    else:
        if bdir != tdir:
            sv_type = "Inversion"
            sv_size = abs(int(bpos) - int(op_pos))
        else:
            if bdir == '+' and int(bpos) <= int(op_pos) or bdir == '-' and int(op_pos) <= int(bpos):
                sv_type = "Deletion"
                sv_size = max(abs(int(bpos) - int(op_pos)) - 1, 0)
            else:
                sv_type = "Tandem Duplication"
                sv_size = max(abs(int(bpos) - int(op_pos)) + 1, 0)


    if sorted([bchr, tchr])[0] == bchr:
        chr1, chr2 = bchr, tchr
    else:
        chr1, chr2 = tchr, bchr

    if sv_type == "Deletion":
        pos1, pos2, dir1, dir2 = min(int(bpos), int(op_pos)), max(int(bpos), int(op_pos)), '+', '-'
    elif sv_type == "Tandem Duplication":
        pos1, pos2, dir1, dir2 = min(int(bpos), int(op_pos)), max(int(bpos), int(op_pos)), '-', '+'
    else:
        if sorted([bchr, tchr])[0] == bchr:
            pos1, pos2, dir1, dir2 = bpos, op_pos, bdir, tdir
        else:
            pos1, pos2, dir1, dir2 = op_pos, bpos, tdir, bdir

    sv_key = ','.join([chr1, str(pos1), dir1, chr2, str(pos2), dir2, str(int(bp_start) - 1)])

    return([sv_key, sv_type, str(sv_size)])


# check whether potential rna contamination or not
def rna_junction_check(chr, start, end, ref_junc_tb, ens_junc_tb, margin = 2):

    junc_flag = False
 
    # check junction annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_junc_tb.fetch(chr, int(start) - 10, int(start) + 10)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            start_diff = int(record[1]) - int(start)
            end_diff = int(record[2]) + 1 - int(end)
            if abs(start_diff - end_diff) <= margin: junc_flag = True

    # check junction annotation for ensGene 
    tabixErrorFlag = 0
    try:
        records = ens_junc_tb.fetch(chr, int(start) - 10, int(start) + 10)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            start_diff = int(record[1]) - int(start)
            end_diff = int(record[2]) + 1 - int(end)
            if abs(start_diff - end_diff) <= margin: junc_flag = True

    return junc_flag




def remove_dup_sv(qname2sv_key):

    tmp_dir = tempfile.mkdtemp()
    qname2dup_flag = {}

    with open(tmp_dir + "/qname2sv_key.tmp.txt", 'w') as hout:
        for qname in qname2sv_key:
            if qname2sv_key[qname] == "---": 
                qname2dup_flag[qname] = "FALSE"
            else:
                # import pdb; pdb.set_trace()
                # print(qname2sv_key[qname])
                tchr1, tend1, tdir1, tchr2, tend2, tdir2, tinseq = qname2sv_key[qname].split(',')
                tstart1, tstart2 = str(int(tend1) - 1), str(int(tend2) - 1)
                print('\t'.join([tchr1, tstart1, tend1, tchr2, tstart2, tend2, qname, tinseq, tdir1, tdir2]), file = hout)

    with open(tmp_dir + "/qname2sv_key.tmp.sorted.txt", 'w') as hout:
        subprocess.call(["sort", "-k1,1", "-k3,3n", "-k4,4", "-k6,6n", tmp_dir + "/qname2sv_key.tmp.txt"], stdout = hout) 

    
    key2info = {}
    with open(tmp_dir + "/qname2sv_key.tmp.sorted.txt", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            del_list = []
            skip_flag = 0

            # if F[2] == "62793856":
            #     import pdb; pdb.set_trace()

            for tkey in key2info:
                # import pdb; pdb.set_trace()
                tchr1, _, tend1, tchr2, _, tend2, _, tinseq, tdir1, tdir2 = key2info[tkey]
                 
                if F[0] != tchr1 or int(F[2]) > int(tend1) + 1000:
                    qname2dup_flag[tkey] = "FALSE"
                    del_list.append(tkey)

                else:
                    if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2 and abs(int(F[2]) - int(tend1)) <= 10 and abs(int(F[5]) - int(tend2)) <= 10:
                        qname2dup_flag[F[6]] = "TRUE"
                        skip_flag = 1

            for tkey in del_list:
                del key2info[tkey]

            if skip_flag == 0:
                key2info[F[6]] = F
   
        # last processing
        for tkey in key2info:
            tchr1, _, tend1, tchr2, _, tend2, _, tinseq, tdir1, tdir2 = key2info[tkey]
            qname2dup_flag[tkey] = "FALSE"
            
    shutil.rmtree(tmp_dir)

    return(qname2dup_flag)




def classify_by_contig_alignment(input_file, output_file, reference_genome, te_seq = None,  simple_repeat = None, remove_rna = None, bwa_option = "-T0 -h300"):

    if remove_rna:
        genome_id, is_grc = check_reference(reference_genome)
        refseq_junc_info = output_file + ".refseq.junc.bed.gz"
        gencode_junc_info = output_file + ".gencode.junc.bed.gz"
        annot_utils.junction.make_junc_info(refseq_junc_info, "refseq", genome_id, is_grc, False)
        annot_utils.junction.make_junc_info(gencode_junc_info, "gencode", genome_id, is_grc, False)
        refseq_junc_tb = pysam.TabixFile(refseq_junc_info)
        gencode_junc_tb = pysam.TabixFile(gencode_junc_info)


    bwa_cmds = ["bwa", "mem"] + bwa_option.split(' ')
    
    key2seq = {}
    hout = open(output_file + ".tmp4.alignment_check.fa", 'w')
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        header2ind = {}
        for (i, cname) in enumerate(header):
            header2ind[cname] = i
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            if F[header2ind["Contig_Post_BP"]] != "---":
                key2seq[key] = F[header2ind["Contig_Post_BP"]]
                print('>' + key + '\n' + F[header2ind["Contig_Post_BP"]], file = hout)

    hout.close()

    FNULL = open(os.devnull, 'w')
    sret = subprocess.call(bwa_cmds + [reference_genome, output_file + ".tmp4.alignment_check.fa",
                           "-o", output_file + ".tmp4.alignment_check.human.sam"], stdout = FNULL, stderr = subprocess.STDOUT)
    if sret != 0:
        print("bwa error, error code: " + str(sret), file = sys.stderr)

    if te_seq is not None:
        sret = subprocess.call(bwa_cmds + [te_seq, output_file + ".tmp4.alignment_check.fa",
                               "-o", output_file + ".tmp4.alignment_check.te.sam"], stdout = FNULL, stderr = subprocess.STDOUT)
        if sret != 0:
            print("bwa error, error code: " + str(sret), file = sys.stderr)

    FNULL.close() 


    samfile = pysam.AlignmentFile(output_file + ".tmp4.alignment_check.human.sam", 'r')

    temp_qname = ''
    temp_read_set = []
    qname2bp_type = {}
    qname2alignment_str = {}
    qname2sv_key = {}
    qname2sv_type = {}
    qname2sv_size = {}

    for read in samfile.fetch():

        if read.query_name != temp_qname:
            if temp_qname != '':
                qname2bp_type[temp_qname] = get_bp_type(temp_read_set)
                bp_pos2alignment = get_alignment_pos(temp_read_set)
                qname2sv_key[temp_qname], qname2sv_type[temp_qname], qname2sv_size[temp_qname] = get_sv_info(temp_qname, bp_pos2alignment[0][0], bp_pos2alignment[0][2]) if len(bp_pos2alignment) >= 1 else ["---", "---", "---"]
                if len(bp_pos2alignment) > 0: 
                    qname2alignment_str[temp_qname] = ';'.join([str(x[0]) + '-' + str(x[1]) + ':' + x[2] for x in bp_pos2alignment])
                else:
                    qname2alignment_str[temp_qname] = "---"

            temp_qname = read.query_name
            temp_read_set = []

        temp_read_set.append(read)

    # last treatment
    if temp_qname != '':

        qname2bp_type[temp_qname] = get_bp_type(temp_read_set)
        bp_pos2alignment = get_alignment_pos(temp_read_set)
        qname2sv_key[temp_qname], qname2sv_type[temp_qname], qname2sv_size[temp_qname] = get_sv_info(temp_qname, bp_pos2alignment[0][0], bp_pos2alignment[0][2]) if len(bp_pos2alignment) >= 1 else ["---", "---", "---"]
        if len(bp_pos2alignment) > 0:
            qname2alignment_str[temp_qname] = ';'.join([str(x[0]) + '-' + str(x[1]) + ':' + x[2] for x in bp_pos2alignment])
        else:
            qname2alignment_str[temp_qname] = "---"

    qname2dup_flag = remove_dup_sv(qname2sv_key)

    samfile.close()

    if te_seq is not None:

        samfile = pysam.AlignmentFile(output_file + ".tmp4.alignment_check.te.sam", 'r')
        temp_qname = ''
        temp_read_seq = []
        qname2alignment_str_te = {}
        
        for read in samfile.fetch():
            if read.query_name != temp_qname:
                if temp_qname != '':
                    bp_pos2alignment = get_alignment_pos(temp_read_set)
                    if len(bp_pos2alignment) > 0:
                        qname2alignment_str_te[temp_qname] = ';'.join([str(x[0]) + '-' + str(x[1]) + ':' + x[2] for x in bp_pos2alignment])
                    else:
                        qname2alignment_str_te[temp_qname] = "---"
                temp_qname = read.query_name
                temp_read_set = []

            temp_read_set.append(read)

        # last treatment
        if temp_qname != '':
            bp_pos2alignment = get_alignment_pos(temp_read_set)
            if len(bp_pos2alignment) > 0:
                qname2alignment_str_te[temp_qname] = ';'.join([str(x[0]) + '-' + str(x[1]) + ':' + x[2] for x in bp_pos2alignment])
            else:
                qname2alignment_str_te[temp_qname] = "---"

    if simple_repeat is not None:
        simple_repeat_tb = pysam.TabixFile(simple_repeat)
 
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        header = header + ["BP_Type", "Human_Alignment", "SV_Key", "SV_Type", "SV_Size", "Is_Dup_SV"]
        if remove_rna: header = header + ["Is_RNA"]
        if te_seq is not None: header = header + ["TE_Alignment"]
        if simple_repeat is not None: header = header + ["Simple_Repeat"]
        print('\t'.join(header), file = hout)

        for (i, cname) in enumerate(header):
            header2ind[cname] = i
        for line in hin:
            F = line.rstrip('\n').split('\t')
            qname = ','.join(F[:4])
            
            bp_type = qname2bp_type[qname] if qname in qname2bp_type else "---"
            alignment_str = qname2alignment_str[qname] if qname in qname2alignment_str else "---"
            sv_key = qname2sv_key[qname] if qname in qname2sv_key else "---"
            sv_type = qname2sv_type[qname] if qname in qname2sv_type else "---"
            sv_size = qname2sv_size[qname] if qname in qname2sv_size else "---"
            dup_flag = qname2dup_flag[qname] if qname in qname2dup_flag else "---"

            print_list = F + [bp_type, alignment_str, sv_key, sv_type, sv_size, str(dup_flag)]

            if remove_rna:
                is_rna = "FALSE"
                if sv_type == "Deletion":
                    tchr1, tpos1, _, _, tpos2, _, _ = sv_key.split(',') 
                    if rna_junction_check(tchr1, tpos1, tpos2, refseq_junc_tb, gencode_junc_tb): is_rna = "TRUE"
                print_list = print_list + [is_rna]

            if te_seq is not None: 
                alignment_str_te = qname2alignment_str_te[qname] if qname in qname2alignment_str_te else "---"
                print_list = print_list + [alignment_str_te]

            
            if simple_repeat is not None:
                simple_repeat_flag = "FALSE"
                tabix_error_flag = 0
                try:
                    records = simple_repeat_tb.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabix_error_flag = 1

                if tabix_error_flag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if int(record[1]) <= int(F[1]) and int(F[1]) <= int(record[2]):
                            simple_repeat_flag = "TRUE"

                print_list = print_list + [simple_repeat_flag]


            print('\t'.join(print_list), file = hout)
            

    subprocess.call(["rm", "-rf", output_file + ".tmp4.alignment_check.fa"])
    subprocess.call(["rm", "-rf", output_file + ".tmp4.alignment_check.human.sam"])
    if remove_rna:
        subprocess.check_call(["rm", "-rf", refseq_junc_info])
        subprocess.check_call(["rm", "-rf", refseq_junc_info + ".tbi"])
        subprocess.check_call(["rm", "-rf", gencode_junc_info])
        subprocess.check_call(["rm", "-rf", gencode_junc_info + ".tbi"])

    if te_seq is not None: subprocess.call(["rm", "-rf", output_file + ".tmp4.alignment_check.te.sam"])
   

 
