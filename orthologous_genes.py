#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  orthologous_genes.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import os
import sys
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from itertools import chain

# import python scripts from the directory of the current script
sdir = f'{os.path.dirname(os.path.realpath(__file__))}/'
sys.path.append('sdir')
import saturation_plot, pangenome_variation


def main(graph_dict, gene_dict, sample_genes, sample_list, output_dir, type_list, cluster_threshold, colour_dict):
    global overlap_dict, ortho_dict
    overlap_dict = {}
    ortho_dict = {}
    # create a dictionary with genes
    parse_graph_genes(graph_dict)
    calc_percentage_overlap()
    # cluster genes separately for each chromosome
    for chromosome, chrom_dict in overlap_dict.items():
        print(f"Calculating orthogroups for {chromosome}")
        cluster_groups = find_best_match(chrom_dict, output_dir, cluster_threshold)
        convert_to_orthogroups(cluster_groups, chromosome, sample_list, sample_genes)
    write_groups(ortho_dict, output_dir, f'orthogroups_{cluster_threshold}')
    # find orthogroup locations for each sample and write them into a file
    location_dict = orthogroup_locations(gene_dict, ortho_dict)
    write_groups(location_dict, output_dir, f'orthogroup_locations_{cluster_threshold}')
    ortho_variation = orthogroup_variation(location_dict, sample_list, type_list)
    # calculate saturation plot
    make_plot(ortho_dict, output_dir, sample_list, type_list, colour_dict, f'orthogroups_{cluster_threshold}')

    return ortho_dict, ortho_variation


def parse_graph_genes(graph_dict):
    for node, node_dict in graph_dict.items():
        gene_list = node_dict['gene_id']
        # parse gene list per chromosome
        chromosome = node_dict['chromosome']
        if chromosome not in overlap_dict:
            overlap_dict[chromosome] = {}
        for (quary_id, quary_region) in gene_list:
            if quary_id not in overlap_dict[chromosome]:
                overlap_dict[chromosome][quary_id] = {}
            # find the overlap between the query and target
            for (target_id, target_region) in gene_list:
                region_overlap = find_overlap(quary_region[0], quary_region[1], target_region[0], target_region[1])
                if region_overlap[1] > region_overlap[0]:
                    overlap_len = region_overlap[1] - region_overlap[0]
                else:
                    overlap_len = 0
                # append dictionary
                if target_id not in overlap_dict[chromosome][quary_id]:
                    overlap_dict[chromosome][quary_id][target_id] = 0
                overlap_dict[chromosome][quary_id][target_id] += overlap_len


def find_overlap(start1, end1, start2, end2):
    # Calculate the start and end positions of the overlap
    overlap_start = max(start1, start2) - start2
    overlap_end = min(end1, end2) - start2

    # Return the start and end position
    return [overlap_start, overlap_end]


def calc_percentage_overlap():
    for chromosome, chrom_dict in overlap_dict.items():
        for query_gene in chrom_dict:
            query_len = overlap_dict[chromosome][query_gene][query_gene]
            for target_gene in overlap_dict[chromosome]:
                if target_gene in overlap_dict[chromosome][query_gene]:
                    target_len = overlap_dict[chromosome][target_gene][target_gene]
                    overlap_len = overlap_dict[chromosome][query_gene][target_gene]
                    perc_overlap = overlap_len / min(query_len, target_len)
                    overlap_dict[chromosome][query_gene][target_gene] = perc_overlap
                    overlap_dict[chromosome][target_gene][query_gene] = perc_overlap


def find_best_match(chrom_dict, output_dir, cluster_threshold):
    df = pd.DataFrame(chrom_dict)
    df.fillna(0, inplace=True)
    # write dataframe to file
    # output_fn = f'{output_dir}/Pe_gene_overlap.tsv'
    # df.to_csv(output_fn, sep='\t')

    # Convert similarity values to distances
    distance_matrix = 1 - df.values
    # Perform hierarchical clustering
    linkage_matrix = linkage(distance_matrix, method='average')
    # Get cluster assignments
    clusters = pd.Series(fcluster(linkage_matrix, cluster_threshold, criterion='distance'),
                         index=df.index, name='Clusters')

    # Create a DataFrame to store the cluster assignments and genes
    result_df = pd.concat([clusters, df.index.to_series()], axis=1)
    # Group genes by cluster
    cluster_groups = result_df.groupby('Clusters')

    return cluster_groups


def convert_to_orthogroups(cluster_groups, chromosome, sample_list, sample_genes):
    # parse clusters
    for cluster_num, cluster_genes in cluster_groups:
        cluster_name = f'{chromosome}_{cluster_num}'
        ortho_dict[cluster_name] = {}
        genes_in_cluster = cluster_genes.iloc[:, 1].tolist()
        # order genes by sample
        for sample in sample_list:
            ortho_dict[cluster_name][sample] = []
            for gene_id in genes_in_cluster:
                gene_sample = sample_genes[gene_id]
                if sample == gene_sample:
                    ortho_dict[cluster_name][sample].append(gene_id)


def write_groups(input_dict, output_dir, name):
    # covert to dataframe
    ortho_df = pd.DataFrame(input_dict).T
    # write dataframe to file
    output_fn = f'{output_dir}/Pe_{name}.tsv'
    ortho_df.to_csv(output_fn, sep='\t')


def orthogroup_locations(gene_dict, ortho_dict):
    location_dict = {}
    for orthogroup, sample_dict in ortho_dict.items():
        location_dict[orthogroup] = {}
        for sample, gene_list in sample_dict.items():
            if gene_list:
                chrom, start, end = gene_dict[gene_list[0]]
                location_dict[orthogroup][sample] = f'{sample}_{chrom}:{start}-{end}'
    return location_dict


def orthogroup_variation(location_dict, sample_list, type_list):
    ortho_variation = {}
    for orthogroup, sample_dict in location_dict.items():
        # extract the samples present in that orthogroup
        ortho_samples = tuple([sample for sample in sample_list if sample in sample_dict])
        # classify orthogroup
        group_variation = pangenome_variation.classify_variation(ortho_samples, sample_list, type_list)
        # append the location and variation for each sample in the group
        for sample in ortho_samples:
            if sample not in ortho_variation:
                ortho_variation[sample] = {}
            chromosome, group_location = sample_dict[sample].split(':')
            if chromosome not in ortho_variation[sample]:
                ortho_variation[sample][chromosome] = {}
            ortho_variation[sample][chromosome][group_location] = group_variation

    return ortho_variation


def make_plot(ortho_dict, output_dir, sample_list, type_list, colour_dict, name):
    counts_dict = count_orthogroups(ortho_dict, sample_list)
    variation_dict = count_variation(counts_dict, sample_list, type_list)
    # write orthogroup variation in a file
    pangenome_variation.write_size(output_dir, variation_dict, name, 'all')
    # plot dataframes to saturation plot
    saturation_plot.plot_saturation(counts_dict, sample_list, output_dir, name, colour_dict)


def count_orthogroups(ortho_dict, sample_list):
    counts_dict = {}
    for orthogroup, sample_dict in ortho_dict.items():
        # extract the samples present in that orthogroup
        ortho_samples = tuple([sample for sample in sample_list if sample_dict[sample]])
        # append count dictionary
        if ortho_samples not in counts_dict:
            counts_dict[ortho_samples] = 0
        counts_dict[ortho_samples] += 1

    return counts_dict


def count_variation(counts_dict, sample_list, type_list):
    variation_dict = {}
    for sample in chain(sample_list, ['Pangenome']):
        variation_dict[sample] = {type_list[0]: 0, type_list[1]: 0, type_list[2]: 0}
        for ortho_samples, count in counts_dict.items():
            # classify orthogroup
            group_variation = pangenome_variation.classify_variation(ortho_samples, sample_list, type_list)
            # append dictionary
            if sample in ortho_samples or sample == 'Pangenome':
                variation_dict[sample][group_variation] += count

    return variation_dict
