#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  saturation_plot.py
#
#  Copyright 2023 Petros Skiadas <p.skiadas@uu.nl>
#

import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import multiprocessing
from functools import partial


def main(graph_dict, sample_list, feature, output_dir, colour_dict):
    counts_dict = count_graph(graph_dict, sample_list, feature)
    plot_saturation(counts_dict, sample_list, output_dir, feature, colour_dict)


def count_graph(graph_dict, sample_list, feature):
    counts_dict = {}
    for node, node_dict in graph_dict.items():
        node_samples = tuple(sample for sample in sample_list if sample in node_dict)
        count = 0
        if isinstance(node_dict[feature], int):
            count = node_dict[feature]
        elif isinstance(node_dict[feature], list):
            count = sum(end - start for start, end in node_dict[feature])
        counts_dict[node_samples] = counts_dict.get(node_samples, 0) + count

    return counts_dict


def calc_variation(sample_combinations, counts_dict, num_processes=None):
    # Determine the number of processes to use
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()

    # Create a multiprocessing Pool with specified number of processes
    pool = multiprocessing.Pool(processes=num_processes)
    # Create a partial function for count_variation with fixed counts_dict
    count_variation_partial = partial(count_variation, counts_dict=counts_dict)
    # Map the count_variation_partial function to each combination in sample_combinations
    results = pool.map(count_variation_partial, sample_combinations)

    # Close the pool to free up resources
    pool.close()
    pool.join()

    # Yield the results
    for result in results:
        yield result


def count_variation(included_samples, counts_dict):
    core_count = 0
    accessory_count = 0
    unique_count = 0
    # find core, accessory, and unique nodes
    for key, count in counts_dict.items():
        key_samples = set(key)
        sub_list_samples = set(included_samples)
        # core if all the samples from the list are found in the key
        if sub_list_samples.issubset(key_samples):
            core_count += count
        # unique if only one sample from the list are found in the key
        elif len(sub_list_samples.intersection(key_samples)) == 1:
            unique_count += count
        # accessory if two or more samples from the list are found in the key
        elif len(sub_list_samples.intersection(key_samples)) >= 2:
            accessory_count += count

    return len(included_samples), core_count, unique_count + accessory_count, accessory_count


def plot_saturation(counts_dict, sample_list, output_dir, feature, colour_dict):
    unique_combinations = generate_combinations(sample_list)
    data = calc_variation(unique_combinations, counts_dict)
    df = pd.DataFrame(data, columns=['Samples', 'Core', 'Unique', 'Accessory'])

    sns.lineplot(data=df, x='Samples', y='Core', label='Core', errorbar='sd', color=colour_dict['core'])
    sns.lineplot(data=df, x='Samples', y='Accessory', label='Accessory', errorbar='sd', color=colour_dict['accessory'])
    sns.lineplot(data=df, x='Samples', y='Unique', label='Accessory + Unique', errorbar='sd', color=colour_dict['unique'])

    # Set x-axis ticks to integers
    plt.xticks(range(int(df['Samples'].min()), int(df['Samples'].max())+1, 2))
    if 'orthogroups' in feature:
        plt.ylabel('Gene orthogroups')
    else:
        y_ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 1e6))
        plt.gca().yaxis.set_major_formatter(y_ticks)
        plt.ylabel('Nucleotides (Mb)')

    plt.xlabel('Included isolates')
    plt.legend()
    plt.savefig(f'{output_dir}/Pe_{feature}_saturation.pdf', format='pdf', bbox_inches='tight', dpi=400)
    plt.close()


def generate_combinations(sample_list):
    # Generate unique combinations of lengths 1 to length of sample list
    unique_combinations = set()
    for r in range(1, len(sample_list) + 1):
        combinations_r = itertools.combinations(sample_list, r)
        for combination in combinations_r:
            unique_combinations.add(tuple(sorted(combination)))

    # Convert the set of tuples back to a list of lists
    unique_combinations = [list(combination) for combination in unique_combinations]
    # append the whole sample list one more time to connect the error bar on the saturation plot
    unique_combinations.append(sample_list)

    return unique_combinations
