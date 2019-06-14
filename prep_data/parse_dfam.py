#! /usr/bin/env python

import sys, gzip

taxon_list = ["Bilateria", "Deuterostomia", "Chordata", "Craniata", "Vertebrata", "Gnathostomata",
            "Teleostomi", "Euteleostomi", "Sarcopterygii", "Dipnotetrapodomorpha", 
            "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria",
            "Euarchontoglires", "Primates", "Haplorrhini", "Simiiformes", "Catarrhini",
            "Hominoidea", "Hominidae", "Homininae", "Homo"]

input_file = sys.argv[1]

cur_tag = ''
cur_id = ''
cur_os = ''
reading_sq = False
cur_seq = ''
cur_nm = ''
with gzip.open(input_file, 'rt') as hin:

    for line in hin:
        line = line.rstrip('\n')
        tag = line[:2]

        if tag == 'ID':
            cur_id = line[2:].split(';')[0].strip()
            if not cur_id.startswith("DF"): 
                print("something is wrong!")
                print(cur_id)
                sys.exit(1)

        if tag == 'NM':
            cur_nm = line[2:].strip()

        if tag == 'OS':
            cur_os = line[2:].strip()

        if tag == 'SQ':
            reading_sq = True
            continue

        if tag == '//':
            if cur_os in taxon_list:
                print('>' + cur_nm)
                print(cur_seq)
            reading_sq = False
            cur_id = ''
            cur_seq = ''
            cur_os = ''

        if reading_sq == True:
            line = line.lstrip() 
            if len(line) >= 65: line = line[:66]
            line = line.replace(' ', '').upper()
            cur_seq = cur_seq + line

