# EBBreak
empirical bayesian framework for detecting somatic breakpoints from sequencing data

## Dependency
### Python (>= 2.7)
Python packages
 - pysam (https://github.com/pysam-developers/pysam, https://pysam.readthedocs.io/en/latest/api.html)
 - annot_utils (https://github.com/friend1ws/annot_utils)
 - numpy
 - scipy
 - statistics

### blat

https://genome.ucsc.edu/cgi-bin/hgBlat,  
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/

### swalign

https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library

### fermi-lite

https://github.com/lh3/fermi-lite

## Install

## Preparation
When using blat -ooc option, prepare ooc file before usage.
```
cd path/to/EBbreak
path/to/faToTwoBit ref/GRCh37.fa ref/GRCh37.2bit
path/to/blat -makeOoc=ref/11.ooc -repMatch=2253 -tileSize=11 ref/GRCh37.2bit /dev/null /dev/null
```

## Command (how to use)
EBbreak consists of five steps.

### parse
```
```

### merge


### filt

### contig


### classify

