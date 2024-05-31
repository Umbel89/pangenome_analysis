#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  orthologous_effectors.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import os
import sys
import glob

# import python scripts from the directory of the current script
sdir = f'{os.path.dirname(os.path.realpath(__file__))}/'
sys.path.append('sdir')
import orthologous_genes


def main(ortho_dict, gene_dict, sample_list, output_dir, type_list, gene_type, input_dir, cluster_threshold, colour_dict):
    effector_list = []
    effector_groups = {}

    parse_effector_list(effector_list, sample_list, input_dir)
    effector_orthogroups(effector_groups, effector_list, ortho_dict)
    # write output of effector orthogroups
    orthologous_genes.write_groups(effector_groups, output_dir, f'{gene_type}_orthogroups_{cluster_threshold}')
    orthologous_genes.make_plot(effector_groups, output_dir, sample_list, type_list, colour_dict,
                                f'{gene_type}_orthogroups_{cluster_threshold}')
    # find orthogroup locations for each sample and write them into a file
    location_dict = orthologous_genes.orthogroup_locations(gene_dict, effector_groups)
    orthologous_genes.write_groups(location_dict, output_dir, f'{gene_type}_orthogroups_locations_{cluster_threshold}')
    ortho_variation = orthologous_genes.orthogroup_variation(location_dict, sample_list, type_list)

    return ortho_variation


def parse_effector_list(effector_list, sample_list, input_dir):
    for sample in sample_list:
        # input tsv files
        gene_tsv = glob.glob(f"{input_dir}/{sample}*.txt")[0]

        with open(gene_tsv) as the_file:
            for line in the_file:
                effector_id = line.split()[0]
                effector_list.append(effector_id)


def effector_orthogroups(effector_groups, effector_list, ortho_dict):
    # if any genes in a cluster are effectors, add the orthogroup
    for cluster_name, cluster_dict in ortho_dict.items():
        for sample, gene_list in cluster_dict.items():
            if any(gene_id in effector_list for gene_id in gene_list):
                effector_groups[cluster_name] = cluster_dict
                break
