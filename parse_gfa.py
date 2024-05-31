#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parse_gfa.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import re
import gzip
import sys
import pandas as pd
from intervaltree import IntervalTree


def main(input_gfa):
    global graph_dict
    global chrom_dict
    global node_dict
    global sample_list
    graph_dict = {}
    chrom_dict = {}
    node_dict = {}
    sample_list = []

    parse_gfa(input_gfa)

    return graph_dict, chrom_dict, sample_list


def parse_gfa(input_gfa):
    with gzip.open(input_gfa, 'rb') as gz_file:
        file_contents = gz_file.read()
        the_file = file_contents.decode('utf-8').split('\n')
        for line in the_file:
            # parse node lines
            if line.startswith('S'):
                node, seq = line.split('\t')[1:3]
                node_dict[node] = len(seq)
            # parse the path for each isolate
            elif line.startswith('W'):
                contig, start, end, graph_path = line.strip().split('\t')[3:]
                sample, chromosome = contig.split('_')
                # populate dictionaries
                if sample not in sample_list:
                    sample_list.append(sample)
                    chrom_dict[sample] = {}
                if chromosome not in chrom_dict[sample]:
                    chrom_dict[sample][chromosome] = IntervalTree()
                # Split the input string at every occurrence of '>' or '<'
                node_path = split_graph_walk(graph_path)
                # append dictionaries
                append_dict(sample, chromosome, node_path, int(start), int(end))


def split_graph_walk(graph_path):
    split_string = re.split('[><]', graph_path)[1:]
    # Replace '>' with '+' and '<' with '-'
    signs = []
    for substring in graph_path:
        if substring == '>':
            signs.append('+')
        elif substring == '<':
            signs.append('-')
    # Combine the signs and numbers into strings
    result = [[sign, item] for sign, item in zip(signs, split_string)]
    return result


def append_dict(sample, chromosome, node_path, start, end):
    region_end = start
    # follow the node path for each sample
    for strand, node in node_path:
        node_length = node_dict[node]
        # calculate end location of each node in the sample chromosome
        region_start = region_end
        region_end += node_length
        # append dictionaries
        if node not in graph_dict:
            graph_dict[node] = {'node_length': node_length, 'gene': [], 'repeat': [], 'gene_id': [],
                                'chromosome': chromosome, 'variation': '', 'samples': ''}
        if sample not in graph_dict[node]:
            graph_dict[node][sample] = []
        graph_dict[node][sample].append([region_start, region_end, strand])
        chrom_dict[sample][chromosome][region_start:region_end] = node

    # as QC, calculate the difference between graph length and chromosome length
    calc_difference(chromosome, end, region_end, sample)


def calc_difference(chromosome, end, region_end, sample):
    difference = end - region_end
    if difference > 0:
        print(f'{sample}_{chromosome}: graph length = {region_end}, chrom length = {end}, difference = {difference}')


def write_output(output_dir, graph_dict):

    df = pd.DataFrame(graph_dict).T
    # write dataframe to file
    output_fn = f'{output_dir}/Pe_graph_nodes.tsv'
    df.to_csv(output_fn, sep='\t')


if __name__ == '__main__':

    input_gfa, output_dir = sys.argv[1:]
    graph_dict, chrom_dict, sample_list = main(input_gfa)
    write_output(output_dir, graph_dict)
