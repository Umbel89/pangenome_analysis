#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parse_gff.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import sys
import re
import glob


def main(sample_list, gene_gff_dir, repeat_gff_dir):
    global gff_dict, gene_dict, sample_genes
    gff_dict = {'gene': {}, 'repeat': {}}
    gene_dict = {}
    sample_genes = {}

    for sample in sample_list:
        # input gff files
        repeat_gff = glob.glob(f"{gene_gff_dir}/{sample}*.gff3")[0]
        gene_gff = glob.glob(f"{repeat_gff_dir}/{sample}*.gff3")[0]

        parse_gff(repeat_gff, sample, 'similarity', 'repeat')
        parse_gff(gene_gff, sample, 'CDS', 'gene')

    return gff_dict, gene_dict, sample_genes


def parse_gff(input_gff, sample, feature, type_key):
    gff_dict[type_key][sample] = {}
    # open input gff and parse only the needed features
    with open(input_gff) as the_file:
        for line in the_file:
            if not line.startswith('#'):
                # feature should be mRNA for genes and similarity for repeats
                gff_feature = line.split()[2]
                gff_chrom = line.split()[0].split('_')[-1]
                if feature == gff_feature and 'Chr18' not in gff_chrom:
                    # find feature location and change the list
                    start, end = (int(x) for x in line.split()[3:5])
                    if type_key == 'gene':
                        trans_id = re.search(r'Parent=(.*?)(;|$)', line).group(1)
                        sample_genes[trans_id] = sample
                        # append the location of the full transcript
                        merge_regions(start, end, gff_chrom, trans_id)
                    else:
                        trans_id = re.search(r'Motif:(.*?)\"', line).group(1)
                    if gff_chrom not in gff_dict[type_key][sample]:
                        gff_dict[type_key][sample][gff_chrom] = []
                    # convert to 0-based coordinates
                    gff_dict[type_key][sample][gff_chrom].append([[start - 1, end], trans_id])


def merge_regions(start, end, gff_chrom, trans_id):
    if trans_id not in gene_dict:
        gene_dict[trans_id] = [gff_chrom, start, end]
    else:
        old_start = gene_dict[trans_id][1]
        gene_dict[trans_id] = [gff_chrom, old_start, end]


if __name__ == '__main__':
    sample_list = sys.argv[1:]

    main(sample_list.split(','))
