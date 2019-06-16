#! /usr/bin/env bash

wget http://dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz

python parse_dfam.py Dfam.embl.gz > dfam.fa

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz

 


