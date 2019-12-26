#! /usr/bin/env python

import gzip
from annot_utils import chr_name

def mod_simple_repeat_info(ucsc_simple_repeat_file, output_file, genome_id, is_grc):

    # create UCSC to GRC chr name corresponding table
    ucsc2grc = {} 
    if is_grc:
        ucsc2grc = chr_name.make_ucsc2grc(genome_id)

    hout = open(output_file, 'w')
    with gzip.open(ucsc_simple_repeat_file, 'rt') as hin:

        for line in hin:

            F = line.rstrip('\n').split('\t')
        
            chr = ucsc2grc[F[1]] if F[1] in ucsc2grc else F[1]
            print(chr + '\t' + '\t'.join(F[2:]), file = hout)

    hout.close()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser("mod_simple_repeat_info")
    parser.add_argument("ucsc_simple_repeat_file", default = None, type = str,
                        help = "the path to simple repeat file downloaded from UCSC site")

    parser.add_argument("output_file", default = None, type = str,
                        help = "the path to output file")

    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                        help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    parser.add_argument("--grc", default = False, action = 'store_true',
                        help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    args = parser.parse_args()

    mod_simple_repeat_info(args.ucsc_simple_repeat_file, args.output_file, args.genome_id, args.grc)


