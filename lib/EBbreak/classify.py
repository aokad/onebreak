#! /usr/bin/env python

from annot_utils.junction import *
from annot_utils.gene import *
from annot_utils.exon import *
import pysam

def classify_canonicalSV(input_file, output_file):

    """
    function for selecting canonical SV from contig results
    """

    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n')

        for line in hin:
            F = line.rstrip('\n').split('\t')
            flag = 0

            # canonical SV definition
            # Left mapping 
            if F[13] != "---":

                # single alignment (human margin equals ---)
                if F[14] == "---":
                    Chr_2 = F[12].split(',')[0]
                    Dir_2 = F[12].split(',')[3]
                    if Dir_2 == '-':
                        Pos_2 = F[12].split(',')[2]
                    if Dir_2 == '+':
                        Pos_2 = F[12].split(',')[1]
                    Inserted_Seq = F[12].split(',')[4]
                    flag = 1

                # best mapping: 0 or 1 mismatch, margin more than 5 bases for both left and right mapping
                elif ((int(F[13]) == 0 or int(F[13]) == 1) and (int(F[14]) >= 5) and (F[17] == "---")) :
                    Chr_2 = F[12].split(',')[0]
                    Dir_2 = F[12].split(',')[3]
                    if Dir_2 == '-':
                        Pos_2 = F[12].split(',')[2]
                    if Dir_2 == '+':
                        Pos_2 = F[12].split(',')[1]
                    Inserted_Seq = F[12].split(',')[4]
                    flag = 1

                elif ((int(F[13]) == 0 or int(F[13]) == 1) and (int(F[14]) >= 5) and (int(F[17]) >= 5)) :
                    Chr_2 = F[12].split(',')[0]
                    Dir_2 = F[12].split(',')[3]
                    if Dir_2 == '-':
                        Pos_2 = F[12].split(',')[2]
                    if Dir_2 == '+':
                        Pos_2 = F[12].split(',')[1]
                    Inserted_Seq = F[12].split(',')[4]
                    flag = 1

            # right mapping
                # single alignment and  no left mapping
            if F[16] != "---":
                if (F[17] == "---" and F[12] == "---"):
                    Chr_2 = F[15].split(',')[0]
                    Dir_2 = F[15].split(',')[3]
                    if Dir_2 == '-':
                        Pos_2 = F[15].split(',')[2]
                    if Dir_2 == '+':
                        Pos_2 = F[15].split(',')[1]
                    Inserted_Seq = F[15].split(',')[4]
                    flag = 1

            if flag == 0: continue

            if F[0] != Chr_2:
                SVtype = "translocation"
            elif F[2] == "+" and Dir_2 == "-": 
                SVtype = "inversion"
            elif F[2] == "-" and Dir_2 == "+":
                SVtype = "inversion"
            elif F[2] == "+" and Dir_2 == "+" and int(F[1]) < int(Pos_2):
                SVtype = "deletion"
            elif F[2] == "-" and Dir_2 == "-" and int(F[1]) > int(Pos_2):
                SVtype = "deletion"
            elif F[2] == "+" and Dir_2 == "+" and int(F[1]) > int(Pos_2): 
                SVtype = "dulplication"
            elif F[2] == "-" and Dir_2 == "-" and int(F[1]) < int(Pos_2):
                SVtype = "dulplication"
            else:
                SVtype = "no SV"

            print >> hout, '\t'.join(F[:3]) + '\t' + Chr_2 + '\t' + Pos_2 + '\t' + Dir_2 + '\t' + Inserted_Seq + '\t' + SVtype + '\t' + '\t'.join(F[4:])

    hin.close()
    hout.close()


def filter_doublecount(input_file, output_file):

    """
    function for filtering double-counted canoical SVs
    """

    hout = open(output_file, 'w')

    breakpoint_dict = {}

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            key1 = F[0] + ',' + F[1] + ',' + F[7]
            key2 = F[3] + ',' + F[4] + ',' + F[7]

            key_margin = []
            flag = 0

            for i in range(-7, 8):
                key_margin.append(F[3] + ',' + str(int(F[4])+i) + ',' + F[7])

            for key in key_margin:
                if key in breakpoint_dict:
                    if breakpoint_dict[key].split(',')[0] == F[0]:
                        if int(breakpoint_dict[key].split(',')[1]) in range(int(F[1]) - 7, int(F[1]) + 8):
                            flag = 1

            if flag == 1:
                continue

            print >> hout, '\t'.join(F)
            breakpoint_dict[key1] = key2

    hin.close()
    hout.close()



def annot_canonicalSV(input_file, output_file, grc, genome_id):

    """
    function for annotating canonical SVs
    """

    make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, grc, False)
    make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, grc, False)

    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")


    hout = open(output_file, 'w')
    header2ind = {}
    print >> hout, '\t'.join(["Chr_1", "Pos_1", "Dir_1", "Chr_2", "Pos_2", "Dir_2", "Inserted_Seq", "Variant_Type", \
                             "Gene_1", "Gene_2", "Exon_1", "Exon_2", "Num_Tumor_Total_Read_Pair", "Num_Tumor_Var_Read_Pair", \
                             "Tumor_VAF", "Num_Control_Ref_Read_Pair", "Num_Control_Var_Read_Pair", "Control_VAF", \
                             "Minus_Log_Fisher_P_value", "Long_Contig_Cap3", "Contig_length",
                             "Junc_Seq_Consistency", "Human_left_Alignment", "Human_left_Mismatch", "Human_left_Margin", "Human_right_Alignment", "Human_right_Mismatch", "Human_right_Margin",
                             "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                             "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                             "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                             "Short_Contig_Cap3", "Contig_length", "Junc_Seq_Consistency", "Human_left_Alignment", "Human_left_Mismatch", "Human_left_argin", "Human_right_Alignment", "Human_right_Mismatch", "Human_right_Margin", \
                             "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin", \
                             "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                             "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                             "Long_Contig_SGA", "Contig_length", "Junc_Seq_Consistency", "Human_left_Alignment", "Human_left_Mismatch", "Human_left_Margin", "Human_right_Alignment", "Human_right_Mismatch", "Human_right_Margin",
                             "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                             "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                             "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin"])

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            ##########
            # check gene annotation for the side 1
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch("chr"+str(F[0]), int(F[1]) - 1, int(F[1]) + 1)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1

            gene1 = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene1.append(record[3])

            if len(gene1) == 0: gene1.append("---")
            gene1 = list(set(gene1))
            ##########

            ##########
            # check gene annotation for the side 2
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch("chr"+str(F[3]), int(F[4]) - 1, int(F[4]))
            except Exception as inst:
                print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorFlag = 1

            gene2 = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene2.append(record[3])

            if len(gene2) == 0: gene2.append("---")
            gene2 = list(set(gene2))
            ##########

            ##########
            # check exon annotation for the side 1
            tabixErrorFlag = 0
            try:
                records = exon_tb.fetch("chr"+str(F[0]), int(F[1]) - 1, int(F[1]))
            except Exception as inst:
                print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorFlag = 1

            exon1 = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    exon1.append(record[3])

            if len(exon1) == 0: exon1.append("---")
            exon1 = list(set(exon1))
            ##########

            ##########
            # check exon annotation for the side 2
            tabixErrorFlag = 0
            try:
                records = exon_tb.fetch("chr"+str(F[3]), int(F[4]) - 1, int(F[4]))
            except Exception as inst:
                print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorFlag = 1
           
            exon2 = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    exon2.append(record[3])

            if len(exon2) == 0: exon2.append("---")
            exon2 = list(set(exon2))
            ##########

            Tumor_VAF = float(F[9]) / float(F[8])
            Tumor_VAF = str(round(Tumor_VAF, 4))

            Control_VAF = float(F[11]) / float(F[10])
            Control_VAF = str(round(Control_VAF, 4))

            print >> hout, '\t'.join(F[:8]) + '\t' + ';'.join(gene1) + '\t' + ';'.join(gene2) + '\t' + ';'.join(exon1) + '\t' + ';'.join(exon2) + '\t' +  \
                           '\t'.join(F[8:10]) + '\t' + str(Tumor_VAF) + '\t' + '\t'.join(F[10:12]) + '\t' + str(Tumor_VAF) + '\t' + '\t'.join(F[12:])

    hin.close()
    hout.close()
    gene_tb.close()
    exon_tb.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])


def classify_non_canonicalSV(input_file, output_file, canonical_sv_file):

    """
    function for selecting non-canonical SV from contig results
    """

    canonical_SV_breakpoint_list = []

    with open(canonical_sv_file, 'r') as hIN:
        #header = hIN.readline().rstrip('\n').split('\t')
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            key1 = F[0] + ',' + F[1] + ',' + F[2]
            if F[5] == "+":
                key2 = F[3] + ',' + F[4] + ',' + "-"
            if F[5] == "-":
                key2 = F[3] + ',' + F[4] + ',' + "+"
            canonical_SV_breakpoint_list.append(key1)
            canonical_SV_breakpoint_list.append(key2)
    hIN.close()


    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        print >> hout, '\t'.join(header)

        for line in hin:
            F = line.rstrip('\n').split('\t')
            # filter canonical SVs
            flag = 0
            for i in range(-5, 6):
                key = F[0] + ',' + str(int(F[1]) + i) + ',' + F[2]
                if key in canonical_SV_breakpoint_list:
                    flag = 1
            if flag == 1:
                continue
            print >> hout, '\t'.join(F)

    hin.close()
    hout.close()



def filter_rna_junction(input_file, output_file, grc, genome_id):

    """
    function for removing potential RNA contamination
    """

    refseq_junc_info = output_file + ".refseq.junc.bed.gz"
    gencode_junc_info = output_file + ".gencode.junc.bed.gz"

    make_junc_info(refseq_junc_info, "refseq", genome_id, grc, False)
    make_junc_info(gencode_junc_info, "gencode", genome_id, grc, False)

    refseq_junc_tb = pysam.TabixFile(refseq_junc_info)
    gencode_junc_tb = pysam.TabixFile(gencode_junc_info)

    hout = open(output_file, 'w')

    with open(input_file, 'r') as hin:
        for line in hin:

            F = line.rstrip('\n').split('\t')
            if F[7] == "deletion":
                if int(F[1]) <= int(F[4]):
                    if junction_check(F[0], F[1], F[4], refseq_junc_tb, gencode_junc_tb): 
                        continue
                else:
                    if junction_check(F[0], F[4], F[1], refseq_junc_tb, gencode_junc_tb): 
                        continue
            print >> hout, '\t'.join(F)


    subprocess.check_call(["rm", "-rf", refseq_junc_info])
    subprocess.check_call(["rm", "-rf", refseq_junc_info + ".tbi"])
    subprocess.check_call(["rm", "-rf", gencode_junc_info])
    subprocess.check_call(["rm", "-rf", gencode_junc_info + ".tbi"])



def junction_check(chr, start, end, ref_junc_tb, ens_junc_tb, margin = 2):

    """
    function for checkiing potential RNA contamination
    """

    junc_flag = False
 
    # check junction annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_junc_tb.fetch("chr"+str(chr), int(start) - 10, int(start) + 10)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            start_diff = int(record[1]) - int(start)
            end_diff = int(record[2]) + 1 - int(end)
            if abs(start_diff - end_diff) <= margin: 
                junc_flag = True

    # check junction annotation for ensGene 
    tabixErrorFlag = 0
    try:
        records = ens_junc_tb.fetch("chr"+str(chr), int(start) - 10, int(start) + 10)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            start_diff = int(record[1]) - int(start)
            end_diff = int(record[2]) + 1 - int(end)
            if abs(start_diff - end_diff) <= margin: 
                junc_flag = True

    return junc_flag



def annotate_break_point(input_file, output_file, grc, genome_id):

    make_gene_info(output_file + ".tmp.refGene.bed.gz", "refseq", genome_id, grc, False)
    make_exon_info(output_file + ".tmp.refExon.bed.gz", "refseq", genome_id, grc, False)

    gene_tb = pysam.TabixFile(output_file + ".tmp.refGene.bed.gz")
    exon_tb = pysam.TabixFile(output_file + ".tmp.refExon.bed.gz")

    hout = open(output_file, 'w')
    header2ind = {}
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        # for (i, cname) in enumerate(header):
        #     header2ind[cname] = i

        print >> hout, '\t'.join(["Chr", "Pos", "Dir", "Junc_Seq", "Gene", "Exon", "Num_Tumor_Total_Read_Pair", "Num_Tumor_Var_Read_Pair",
                                 "Num_Control_Ref_Read_Pair", "Num_Control_Var_Read_Pair", "Minus_Log_Fisher_P_value", "Long_Contig_Cap3", "Contig_length",
                                 "Junc_Seq_Consistency", "Human_left_Alignment", "Human_left_Mismatch", "Human_left_Margin", "Human_right_Alignment", "Human_right_Mismatch", "Human_right_Margin",
                                 "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                                 "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                 "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                                 "Short_Contig_Cap3", "Contig_length", "Junc_Seq_Consistency", "Human_left_Alignment", "Human_left_Mismatch", "Human_left_argin", "Human_right_Alignment", "Human_right_Mismatch", "Human_right_Margin", \
                                 "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin", \
                                 "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                 "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin",
                                 "Long_Contig_SGA", "Contig_length", "Junc_Seq_Consistency", "Human_left_Alignment", "Human_left_Mismatch", "Human_left_Margin", "Human_right_Alignment", "Human_right_Mismatch", "Human_right_Margin",
                                 "Virus_Alignment", "Virus_Mismatch", "Virus_Margin", "Repeat_Alignment", "Repeat_Mismatch", "Repeat_Margin",
                                 "Mitochondria_Alignment", "Mitochondria_Mismatch", "Mitochondria_Margin",
                                 "Adapter_Alignment", "Adapter_Mismatch", "Adapter_Margin"])

        for line in hin:
            F = line.rstrip('\n').split('\t')

            ##########
            # check gene annotation
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch("chr"+str(F[0]), int(F[1]) - 1, int(F[1]) + 1)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1

            gene = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene.append(record[3])

            gene = list(set(gene))
            if len(gene) == 0: gene.append("---")  

            ##########
            # check gene annotation
            tabixErrorFlag = 0
            try:
                records = exon_tb.fetch("chr"+str(F[0]), int(F[1]) - 1, int(F[1]) + 1)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag = 1
                
            exon = [];
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    exon.append(record[3])
                    
            exon = list(set(exon))
            if len(exon) == 0: exon.append("---")

            print >> hout, '\t'.join(F[:4]) + '\t' + \
                           ','.join(gene) + '\t' + ';'.join(exon) + '\t' + '\t'.join(F[4:])

    hin.close()
    hout.close()
    gene_tb.close()
    exon_tb.close()

    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refGene.bed.gz.tbi"])
    subprocess.call(["rm", "-rf", output_file + ".tmp.refExon.bed.gz.tbi"])
    