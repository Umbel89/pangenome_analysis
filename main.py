#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  main.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import os
import sys
import argparse

# import python scripts from the directory of the current script
sdir = f'{os.path.dirname(os.path.realpath(__file__))}/'
sys.path.append('sdir')
import parse_gfa, parse_gff, saturation_plot, pangenome_variation, orthologous_genes, \
       orthologous_effectors, plot_orthogroups


def integrate_gff():
    for ([chrom_start, chrom_end], trans_id) in region_list:
        for interval in chrom_dict[sample][chromosome][chrom_start:chrom_end]:
            node = interval.data
            # find overlap of the feature with the node region
            region_overlap = orthologous_genes.find_overlap(chrom_start, chrom_end, interval.begin, interval.end)
            # merge overlap with existing overlaps in this node
            add_region_to_node(node, feature, region_overlap, trans_id)


def add_region_to_node(node, feature, region, trans_id):
    region_list = graph_dict[node][feature]
    # Check for overlap and merge regions if necessary
    merged = False
    for i in range(len(region_list)):
        if region_list[i][1] >= region[0] and region_list[i][0] <= region[1]:
            region_list[i] = [min(region_list[i][0], region[0]), max(region_list[i][1], region[1])]
            merged = True
            break

    # If no overlap, add the region separately
    if not merged:
        # Check if the new region extends beyond the previously merged regions
        non_overlapping = True
        for i in range(len(region_list)):
            if region[1] < region_list[i][0] or region[0] > region_list[i][1]:
                continue
            else:
                non_overlapping = False
                break

        # Add the region separately if it does not overlap
        if non_overlapping:
            region_list.append(region)

    # Process the regions until no further overlaps can be merged
    merge_occurred = True
    while merge_occurred:
        merge_occurred = False
        for i in range(len(region_list) - 1):
            if region_list[i][1] >= region_list[i + 1][0]:
                region_list[i] = [min(region_list[i][0], region_list[i + 1][0]),
                                  max(region_list[i][1], region_list[i + 1][1])]
                del region_list[i + 1]
                merge_occurred = True
                break

    # Update the graph dictionary with the final merged regions for the node and append the gene ID
    graph_dict[node][feature] = region_list
    if feature == 'gene' and trans_id not in graph_dict[node]['gene_id']:
        graph_dict[node]['gene_id'].append([trans_id, region])


def input_file(string):
    """Check if file exists."""
    if string != "" and string is not None:
        string = os.path.abspath(string)
        if not os.path.isfile(string):
            print("The input is not a file.")
            raise IOError("not a file [%s]" % string)

    return string


def make_directory(directory):
    if not os.path.exists(os.path.dirname(directory+'/')):
        os.makedirs(os.path.dirname(directory+'/'))

    return directory


def input_options():
    # default options
    def_input_chrom = 'all'
    def_cluster_threshold = 0.6
    # input parser
    parser = argparse.ArgumentParser(description="Parse and annotate a pangenome graph, and plot its variation.")
    parser.add_argument('-i', '--input_gfa', metavar="FILE", type=input_file, required=True,
                        help="File location of pangenome graph gfa.gz file.")
    parser.add_argument('-g', '--gene_gff_dir', metavar="STR", type=str, required=True,
                        help="Directory with the gene gff files for all samples in the graph. dir/[sample]*.gff3")
    parser.add_argument('-r', '--repeat_gff_dir', metavar="STR", type=str, required=True,
                        help="Directory with the repeat gff files for all samples in the graph. dir/[sample]*.gff3")
    parser.add_argument('-e', '--effector_dir', metavar="STR", type=str, required=True,
                        help="Directory with txt files of a subgroup of genes, one gene_id per line, for all samples in"
                             " the graph. dir/[sample]*.txt")
    parser.add_argument('-o', '--output_dir', metavar="STR", type=make_directory, required=True,
                        help="Directory where output will be written.")
    # optional arguments
    parser.add_argument('-c', '--input_chrom', default=def_input_chrom, type=str, metavar="STR",
                        help=f"Specify a chromosome to be parsed. [default={def_input_chrom}]")
    parser.add_argument('-t', '--cluster_threshold', default=def_cluster_threshold, type=float, metavar="FLOAT",
                        help=f"Define orthogroups of genes that their distance is bellow this threshold."
                             f" [default={def_cluster_threshold}]")

    return parser.parse_args()


if __name__ == '__main__':

    args = input_options()
    type_list = ['core', 'accessory', 'unique']
    colour_dict = {type_list[0]: (64/255, 161/255, 201/255),
                   type_list[1]: (146/255, 208/255, 80/255),
                   type_list[2]: (254/255, 192/255, 15/255)}

    # parse input files
    print("Parsing input pangenome graph")
    graph_dict, chrom_dict, sample_list = parse_gfa.main(args.input_gfa)
    print("Annotating pangenome graph")
    gff_dict, gene_dict, sample_genes = parse_gff.main(sample_list, args.gene_gff_dir, args.repeat_gff_dir)
    # overlap annotations with the graph
    for feature in ['gene', 'repeat']:
        for sample, feature_dict in gff_dict[feature].items():
            for chromosome, region_list in feature_dict.items():
                if chromosome == args.input_chrom or args.input_chrom == 'all':
                    integrate_gff()

    # calculate variation
    print("Calculating and plotting pangenome variation")
    graph_dict = pangenome_variation.main(graph_dict, sample_list, args.output_dir, args.input_chrom, type_list)
    # write graph dictionary to a file for manual curation
    # parse_gfa.write_output(args.output_dir, graph_dict)

    # produce saturation plots on the nucleotide level
    for feature in ['repeat', 'node_length']:
        saturation_plot.main(graph_dict, sample_list, feature, args.output_dir, colour_dict)

    # calculate the best matches between genes
    print("Calculating gene orthogroups and plotting gene variation")
    ortho_dict, gene_variation = orthologous_genes.main(graph_dict, gene_dict, sample_genes, sample_list,
                                                        args.output_dir, type_list, args.cluster_threshold, colour_dict)
    effector_variation = orthologous_effectors.main(ortho_dict, gene_dict, sample_list, args.output_dir, type_list,
                                                    'effector', args.effector_dir, args.cluster_threshold, colour_dict)
    # plot orthogroup variation on the chromosomes
    if args.input_chrom == 'all':
        plot_orthogroups.main(gene_variation, sample_list, args.output_dir, colour_dict,
                              f"gene_variation_{args.cluster_threshold}")
        plot_orthogroups.main(effector_variation, sample_list, args.output_dir, colour_dict,
                              f"effector_variation_{args.cluster_threshold}")
