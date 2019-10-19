#! /usr/bin/env python

import sys


def check_reference(reference_genome):

    hg19_dict = {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "X": 155270560,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "20": 63025520,
        "Y": 59373566,
        "19": 59128983,
        "22": 51304566,
        "21": 48129895
    }

    hg38_dict = {
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
        "6": 170805979,
        "7": 159345973,
        "X": 156040895,
        "8": 145138636,
        "9": 138394717,
        "11": 135086622,
        "10": 133797422,
        "12": 133275309,
        "13": 114364328,
        "14": 107043718,
        "15": 101991189,
        "16": 90338345,
        "17": 83257441,
        "18": 80373285,
        "20": 64444167,
        "19": 58617616,
        "Y": 57227415,
        "22": 50818468,
        "21": 46709983
    }

    target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

    ref2length = {}
    tmp_ref, tmp_length = "", 0
    with open(reference_genome) as hin:
        for line in hin:
            if line.startswith('>'):
                if tmp_ref != "":
                    ref2length[tmp_ref] = tmp_length
                tmp_ref = line.replace('>', '').rstrip('\n').split(' ')[0]
                tmp_length = 0
            else:
                tmp_length = tmp_length + len(line.rstrip('\n'))
                
    if tmp_ref != "":
        ref2length[tmp_ref] = tmp_length

    grc_count, nongrc_count = 0, 0
    hg19_count, hg38_count = 0, 0

    for rname in ref2length:
        seq_length = ref2length[rname]

        if rname in target_rnames:
            grc_count = grc_count + 1
            # seq_length = f.get_reference_length(rname)
        elif rname in ["chr" + x for x in target_rnames]:
            nongrc_count = nongrc_count + 1
            # seq_length = f.get_reference_length(rname)
            rname = rname.replace("chr", '')
        else:
            continue

        if seq_length == hg19_dict[rname]:
            hg19_count = hg19_count + 1
        elif seq_length == hg38_dict[rname]:
            hg38_count = hg38_count + 1
        else:
            print("The length of %s is not included in either hg19 or hg38 in the BAM file." % rname, file = sys.stderr)


    if grc_count == 0 and nongrc_count == 0:
        print("No human chromosome is found in the BAM file.", file = sys.stderr)
        sys.exit(1)

    if hg19_count == 0 and hg38_count == 0:
        print("Neither hg19 nor hg38 chromosome is included in the BAM file.", file = sys.stderr)
        sys.exit(1)
 
    if grc_count > 0 and nongrc_count > 0:
        print("Both UCSC (starting with 'chr') and GRC chromosome nomenclature is mixed in the BAM file.", file = sys.stderr)
        sys.exit(1)

    if hg19_count > 0 and hg38_count > 0:
        print("Both hg19 and hg38 chromosome is included in the BAM file.", file = sys.stderr)
        sys.exit(1)

    genome_id = "hg19" if hg19_count > hg38_count else "hg38"
    is_grc = True if grc_count > nongrc_count else False

    return([genome_id, is_grc])



if __name__ == "__main__":

    import sys
    print(check_reference(sys.argv[1]))

