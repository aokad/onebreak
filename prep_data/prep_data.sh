#! /usr/bin/env bash

set -x 
set -e


wget http://dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz

python parse_dfam.py Dfam.embl.gz > dfam.fa

# wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz

python3 mod_simple_repeat_info.py simpleRepeat.txt.gz simple_repeat.bed.gz.unsorted.tmp --grc

sort -k1,1 -k2,2n -k3,3n simple_repeat.bed.gz.unsorted.tmp > simple_repeat.bed.gz.sorted.tmp  

bgzip -f -c simple_repeat.bed.gz.sorted.tmp > simple_repeat.bed.gz

tabix -p bed simple_repeat.bed.gz

rm -rf simple_repeat.bed.gz.unsorted.tmp
rm -rf simple_repeat.bed.gz.sorted.tmp  


