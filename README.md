# Pangenome graph analysis
Parse a pangenome graph and overlap given annotations of genes and TEs for every sample in the graph. For all genes, define orthogroups based on their synteny on the pangenome graph. Plot the variation for the whole pangenome and for repeats at the nucleotide level, and for all genes and a given subset of genes at the orthogroup level.

![Pan_orthogroups-01](https://github.com/Umbel89/pangenome_analysis/assets/94618033/d91e36cc-cb86-4dbd-b886-091e4036514c)

## Installation
```bash
git clone https://github.com/Umbel89/pangenome_analysis.git
```

## Usage
```
python main.py --help

usage: main.py [-h] -i FILE -g STR -r STR -e STR -o STR [-c STR] [-t FLOAT]

Parse and annotate a pangenome graph, and plot its variation.

  -h, --help            show this help message and exit

required arguments:
  -i FILE, --input_gfa FILE
                        File location of pangenome graph gfa.gz file.
  -g STR, --gene_gff_dir STR
                        Directory with the gene gff files for all samples in the graph. dir/[sample]*.gff3
  -r STR, --repeat_gff_dir STR
                        Directory with the repeat gff files for all samples in the graph. dir/[sample]*.gff3
  -e STR, --effector_dir STR
                        Directory with txt files of a subgroup of genes, one gene_id per line, for all samples in the
                        graph. dir/[sample]*.txt
  -o STR, --output_dir STR
                        Directory where output will be written.

optional arguments:
  -c STR, --input_chrom STR
                        Specify a chromosome to be parsed. [default=all]
  -t FLOAT, --cluster_threshold FLOAT
                        Define orthogroups of genes that their distance is bellow this threshold. [default=0.6]

```

## Cite
[Skiadas, P. *et al*. (2024), Unpublished: Pangenome graph analysis reveals extensive effector copy-number variation in spinach downy mildew](https://www.biorxiv.org/content/10.1101/2024.05.30.596583v1)
