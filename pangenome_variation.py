#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pangenome_variation.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import pandas as pd
from itertools import chain


def main(graph_dict, sample_list, output_dir, chromosome, type_list):
    global feature_list, size_list, desired_order
    feature_list = ['node_length', 'gene', 'repeat']
    size_list = ['SNP', 'short', 'long']

    populate_dictionary(sample_list, type_list)
    append_dictionary(graph_dict, sample_list, type_list)
    # write dictionaries to file
    write_variation(output_dir, chromosome)
    write_size(output_dir, size_dict, 'size', chromosome)
    write_size(output_dir, size_count_dict, 'size_count', chromosome)

    return graph_dict


def populate_dictionary(sample_list, type_list):
    global variation_dict, size_dict, size_count_dict
    variation_dict = {}
    size_dict = {}
    size_count_dict = {}
    # calculate variation per sample
    for sample in chain(sample_list, ['Pangenome']):
        variation_dict[sample] = {}
        size_dict[sample] = {}
        size_count_dict[sample] = {}
        # calculate per feature
        for feature in feature_list:
            # calculate per variation type
            for variation_type in type_list:
                variation_dict[sample][f'{feature}_{variation_type}'] = 0
                # calculate per node length
                for size in size_list:
                    size_dict[sample][f'{feature}_{variation_type}_{size}'] = 0
                    size_count_dict[sample][f'{feature}_{variation_type}_{size}'] = 0


def append_dictionary(graph_dict, sample_list, type_list):
    # for every node find the variation
    for node, node_dict in graph_dict.items():
        variation_type, region_list = calc_variation(sample_list, node_dict, type_list)
        # for every feature calculate the length
        for feature in feature_list:
            count = calculate_length(feature, node_dict)
            # append for every sample that is present in this node
            for sample in chain(region_list, ['Pangenome']):
                variation_dict[sample][f'{feature}_{variation_type}'] += count
                # classify node size
                size = classify_size(count)
                if size:
                    size_dict[sample][f'{feature}_{variation_type}_{size}'] += 1
                    size_count_dict[sample][f'{feature}_{variation_type}_{size}'] += count


def classify_size(count):
    if count == 1:
        size = size_list[0]
    elif count > 50:
        size = size_list[2]
    elif 1 < count < 50:
        size = size_list[1]
    else:
        size = None

    return size


def calculate_length(feature, node_dict):
    # extract length of the feature and count
    count = 0
    # sometimes the value is an integer and others it's a list of regions
    if isinstance(node_dict[feature], int):
        count += node_dict[feature]
    elif isinstance(node_dict[feature], list):
        for [start, end] in node_dict[feature]:
            count += end - start

    return count


def calc_variation(sample_list, node_dict, type_list):
    # create a list of all the regions in the dictionary
    region_list = [sample for sample in sample_list if sample in node_dict]
    variation_type = classify_variation(region_list, sample_list, type_list)
    # append dictionary
    node_dict['variation'] = variation_type
    node_dict['samples'] = ','.join(region_list)

    return variation_type, region_list


def classify_variation(region_list, sample_list, type_list):
    if len(region_list) == 1:
        variation_type = type_list[2]
    elif len(region_list) == len(sample_list):
        variation_type = type_list[0]
    else:
        variation_type = type_list[1]
    return variation_type


def write_variation(output_dir, chromosome):
    df = pd.DataFrame(variation_dict).T
    # calculate non_coding regions for each category
    df['non_coding_core'] = df['node_length_core'] - df['gene_core'] - df['repeat_core']
    df['non_coding_accessory'] = df['node_length_accessory'] - df['gene_accessory'] - df['repeat_accessory']
    df['non_coding_unique'] = df['node_length_unique'] - df['gene_unique'] - df['repeat_unique']
    # write dataframe to file
    output_fn = f'{output_dir}/Pe_{chromosome}_pangenome_variation.tsv'
    df.to_csv(output_fn, sep='\t')


def write_size(output_dir, input_dict, output_name, chromosome):
    df = pd.DataFrame(input_dict).T
    # write dataframe to file
    output_fn = f'{output_dir}/Pe_{chromosome}_pangenome_{output_name}.tsv'
    df.to_csv(output_fn, sep='\t')
