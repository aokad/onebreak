import os
import subprocess

def make_fasta_file(input_file, output_file, seq_len_file):
    hout = open(output_file, 'w')
    hout2 = open(seq_len_file, 'w')
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        header2ind = {}
        for (i, cname) in enumerate(header):
            header2ind[cname] = i
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join(F[:4])
            if F[header2ind["Contig_Post_BP"]] != "---":
                print('>' + key + '\n' + F[header2ind["Contig_Post_BP"]], file = hout)
                print(key + '\t' + str(len(F[header2ind["Contig_Post_BP"]])), file = hout2)
    hout2.close()
    hout.close()

def proc_rmsk_info(total_len, rmsk_info):
    alu_len, alu_count = 0, 0
    L1_len, L1_count, L1_dir = 0, 0, []
    SVA_len, SVA_count = 0, 0
    repeat_len, repeat_count = 0, 0
    polyAT_len = 0

    repeat_info = []
    for i in range(len(rmsk_info)):

        F = rmsk_info[i]

        if F[8] == 'C': F[8] = '-'
        if F[8] == '+':
            repeat_key = ','.join([F[5], F[6], F[8], F[9], F[10], F[11], F[12]])
        else:
            repeat_key = ','.join([F[5], F[6], F[8], F[9], F[10], F[12], F[13]])
        repeat_info.append(repeat_key)

        if F[10] == "LINE/L1":
            L1_len = L1_len + int(F[6]) - int(F[5]) + 1
            L1_count = L1_count + 1
            L1_dir.append(F[8])
        if F[10] == "SINE/Alu":
            alu_len = alu_len + int(F[6]) - int(F[5]) + 1
            alu_count = alu_count + 1
        if F[10] == "Retroposon/SVA":
            SVA_len = SVA_len + int(F[6]) - int(F[5]) + 1
            SVA_count = SVA_count + 1
        if F[9] == "(T)n" and i == 0:
            polyAT_len = polyAT_len + int(F[6])   
        if F[9] == "(A)n" and i == len(rmsk_info) - 1:
            polyAT_len = polyAT_len + total_len - int(F[5]) + 1


    repeat_class = "None"
    L1_ratio, Alu_ratio, SVA_ratio = 0.0, 0.0, 0.0
    if float(total_len) - float(polyAT_len) > 0:
        L1_ratio = min(1.0, float(L1_len) / (float(total_len) - float(polyAT_len)))
        Alu_ratio = min(1.0, float(alu_len) / (float(total_len) - float(polyAT_len)))
        SVA_ratio = min(1.0, float(SVA_len) / (float(total_len) - float(polyAT_len)))


    if L1_ratio >= 0.8:
        if L1_count == 1:
            repeat_class = "Simple_LINE1"
        elif L1_count == 2 and L1_dir[0] != L1_dir[1]:
            repeat_class = "Inverted_LINE1"
        else:
            repeat_class = "Other_LINE1"

    if Alu_ratio >= 0.8:
        repeat_class = "Alu"

    if SVA_ratio >= 0.8:
        repeat_class = "SVA"

    return([repeat_class, str(round(L1_ratio, 4)), str(round(Alu_ratio, 4)), str(round(SVA_ratio, 4)), repeat_info])

def summarize_rmsk(input_file, seq_len_file, output_file):

    seq_len = {}
    with open(seq_len_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split()
            if len(F) <= 1: continue
            seq_len[F[0]] = int(F[1])

    temp_key = ''
    temp_rmsk_info = []

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        #print('\t'.join(['Key', 'Repeat_Class', 'L1_Ratio', 'Alu_Ratio', 'SVA_Ratio', 'Repeat_Info']), file = hout)
        for line in hin:
            F = line.rstrip('\n').split()
            if len(F) <= 1 or F[0] in ["SW", "score"]: continue
            if F[0] == "There": continue

            if temp_key != F[4]:
                if temp_key != '':

                    repeat_class, L1_ratio, Alu_ratio, SVA_ratio, repeat_info = proc_rmsk_info(seq_len[temp_key], temp_rmsk_info)
                    print(temp_key + '\t' + repeat_class + '\t' + L1_ratio + '\t' + Alu_ratio + '\t' + SVA_ratio + '\t' + ';'.join(repeat_info), file = hout)

                temp_key = F[4]
                temp_rmsk_info = []
            
            temp_rmsk_info.append(F)

        if temp_key != '':
            repeat_class, L1_ratio, Alu_ratio, SVA_ratio, repeat_info = proc_rmsk_info(seq_len[temp_key], temp_rmsk_info)
            print(temp_key + '\t' + repeat_class + '\t' + L1_ratio + '\t' + Alu_ratio + '\t' + SVA_ratio + '\t' + ';'.join(repeat_info), file = hout)

def insert_classify_main(sv_list_file, output_file, debug = False):

    make_fasta_file(sv_list_file, output_file + ".tmp.fasta", output_file + ".tmp.seq_len.txt")

    ##########
    # repeat masker
    tmp_dir = os.path.dirname(output_file) + "/RepeatMasker.tmp"
    os.makedirs(tmp_dir, exist_ok=True)
    subprocess.check_call(["RepeatMasker", "-species", "human", output_file + ".tmp.fasta", "-dir", tmp_dir])

    summarize_rmsk(tmp_dir + '/' + os.path.basename(output_file + ".tmp.fasta") + ".out", output_file + ".tmp.seq_len.txt", output_file + ".tmp.rmsk.txt")

    rmsk = {}
    with open(output_file + ".tmp.rmsk.txt", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split()
            if len(F) <= 1: continue
            rmsk[F[0]] = F[1:]
    
    hout = open(output_file, 'w')
    with open(sv_list_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        header = header + ['Repeat_Class', 'L1_Ratio', 'Alu_Ratio', 'SVA_Ratio', 'Repeat_Info']
        print('\t'.join(header), file = hout)

        for line in hin:
            F = line.rstrip('\n').split('\t')
            qname = ','.join(F[:4])
            if qname in rmsk:
                print_list = F + rmsk[qname]
            else:
                print_list = F + ["---"] * 5
            print('\t'.join(print_list), file = hout)

    if debug == False:
        subprocess.call(["rm", "-rf", output_file + ".tmp.fasta"])
        subprocess.call(["rm", "-rf", output_file + ".tmp.seq_len.txt"])
        subprocess.call(["rm", "-rf", output_file + ".tmp.rmsk.txt"])
        subprocess.call(["rm", "-rf", tmp_dir])

if __name__ == '__main__':
    #import sys
    #sv_list_file = sys.argv[1]
    #output_file = sys.argv[2]

    sv_list_file = "/home/aiokada/sandbox/onebreak/output/TCGA-ESCA_0.1.0b7/TCGA-2H-A9GG-01A-11R-A37I-31/TCGA-2H-A9GG-01A-11R-A37I-31.onebreak-contig.txt"
    output_file = "/home/aiokada/sandbox/onebreak/output/TCGA-ESCA_0.1.0b7/TCGA-2H-A9GG-01A-11R-A37I-31/TCGA-2H-A9GG-01A-11R-A37I-31.onebreak-mei.txt"
    insert_classify_main(sv_list_file, output_file, debug = True)

    #for sample in ["TCGA-2H-A9GG-01A-11R-A37I-31", "TCGA-2H-A9GN-01A-11R-A37I-31", "TCGA-IC-A6RE-01A-11R-A336-31", "TCGA-IG-A3QL-01A-11R-A24K-31", "TCGA-IG-A50L-01A-11R-A260-31", "TCGA-L5-A4OS-01A-11R-A28J-31", "TCGA-L5-A8NF-01A-11R-A37I-31", "TCGA-LN-A4A1-01A-21R-A260-31", "TCGA-R6-A8W5-01B-11R-A37I-31", "TCGA-V5-AASX-11A-11R-A38D-31"]:
    #
    #    sv_list_file = "/home/aiokada/sandbox/onebreak/output/TCGA-ESCA_0.1.0b7/{sample}/{sample}.onebreak-contig.txt".format(sample = sample)
    #    output_file = "/home/aiokada/sandbox/onebreak/output/TCGA-ESCA_0.1.0b7/{sample}/{sample}.onebreak-classify-RepeatMasker.txt".format(sample = sample)
    #    insert_classify_main(sv_list_file, output_file)

