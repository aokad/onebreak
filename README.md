# onebreak
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

### bwa

https://github.com/lh3/bwa

### fermi-lite

https://github.com/lh3/fermi-lite

### swalign

https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library

## Command (how to use)
onebreak consists of five steps.

### workflow

![image](https://user-images.githubusercontent.com/13672949/161904048-be8bf771-73c5-452a-849d-5696791b2b0f.png)

### parse
```
onebreak parse \
    ${bam_file} \
    ${parse_file} \
    --reference_genome ${reference}
```
The extension of the ${bam_file} is either .cram or .bam.  
The extension of the ${parse_file} is .txt.gz.  

### merge
```
ls /path/to/onebreak-parse/*.txt.gz > ./merge_input_list.txt
onebreak merge_control \
    ./merge_input_list.txt \
    ${merged_control}
```
The extension of the ${merged_control} is .txt.gz.

### filt
```
onebreak filt \
  --merged_control_file ${merged_control} \
  ${parse_file} \
  ${bam_file} \
  ${filt_file} \
  ${reference}
```
The extension of the ${filt_file} is .txt.gz.

### contig
```
onebreak contig \
    ${filt_file} \
    ${bam_file} \
    ${contig_file} \
    ${reference}
```
The extension of the ${contig_file} is .txt.gz.

### classify
```
onebreak classify \
    ${contig_file} \
    ${classify_file} \
    ${reference}
```
The extension of the ${classify_file} is .txt.
