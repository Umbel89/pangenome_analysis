#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_orthogroups.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import re
import seaborn as sns
import matplotlib.pyplot as plt


def main (ortho_variation, sample_list, output_dir, colour_dict, output_name):

    for sample in sample_list:
        genome_fai = f"/net/virus/linuxhome/michael-group/petros/Pe_genome_assemblies/" \
                     f"{sample}_nanopore/final/{sample}_assembly.fasta.fai"
        chromosome_sizes = parse_chromosome_sizes(genome_fai)
        draw_chromosomes_with_genes(sample, chromosome_sizes, ortho_variation[sample], output_dir, colour_dict, output_name)


def parse_chromosome_sizes(fai_file):
    chromosome_sizes = {}
    with open(fai_file) as fai:
        for line in fai:
            parts = line.strip().split("\t")
            chromosome_name, size = parts[0], int(parts[1])
            if 'Chr18' not in chromosome_name:
                chromosome_sizes[chromosome_name] = size
    # Sort the dictionary based on the keys with natural sorting
    sorted_dict = {k: chromosome_sizes[k] for k in sorted(chromosome_sizes, key=natural_sort_key)}
    return sorted_dict


def natural_sort_key(dict_key):
    parts = re.split('(\d+)', dict_key)  # Split the key into non-digit and digit parts
    return [int(part) if part.isdigit() else part for part in parts]


def draw_chromosomes_with_genes(sample, chromosome_sizes, gene_data, output_dir, colour_dict, output_name):
    # Create a Seaborn plot
    sns.set(style="whitegrid", rc={"axes.facecolor": (0, 0, 0, 0)})
    plt.figure(figsize=(8, 6))

    # Create a Seaborn barplot to represent chromosomes
    ax = sns.barplot(x=[size for chromosome, size in chromosome_sizes.items()],
                     y=chromosome_sizes.keys(), color="white", edgecolor="grey", linewidth=0.3)

    # Customize gene lines (you can adjust the heights for positioning)
    for chromosome in chromosome_sizes:
        if chromosome in gene_data:
            for gene_location, variation in gene_data[chromosome].items():
                start, end = map(int, gene_location.split('-'))
                gene_y = list(chromosome_sizes.keys()).index(chromosome) # Position the gene line at the chromosome's height
                plt.hlines(y=gene_y, xmin=start, xmax=end, color=colour_dict[variation], lw=12)

    ax.set(xlabel="Chromosome Position (Mb)", title=(output_name.title()))
    # Save figure
    plt.savefig(f'{output_dir}/{sample}_{output_name}.pdf', format='pdf', bbox_inches='tight', dpi=400)
    plt.close()
