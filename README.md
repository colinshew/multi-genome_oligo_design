# multi-genome_oligo_design
Identify homologous 200mers within and between genomes

## Description:

This Snakefile takes a BED input and finds homologous 200mer sequences from these regions (tiled at 2x coverage) within and between genomes. These sequences are designed for use in massively parellel reporter assays.

In brief, intervals are merged, divided into overlapping 200-base windows, and the underlying sequence is extracted. Each "tile" is assigned a strand based on the nearest transcribed feature and aligned to each provided genome with BLAT. Alignments are filtered in two steps: a first pass keeps oligos with fewer N good alignments, and a second pass can be tuned to be more stringent (this rescues reads from duplicated regions that have multiple "good" alignments but may still generate spurious hits; the desired alignment percent identity, length, and number of hits can be set for each genome). Finally, only oligos with alignments passing the filter in all genomes are retained. These sequences, as well as a set of scrambled controls, are given universal priming sequences and screened for restriction enzyme cut sites.

## Dependencies

Snakefile developed with the following:
- BLAT (0.35)
- Bedtools (2.25.0) 
- Python (3.6.7)

The the following Python libraries are also used:
- SeqIO (1.72)
- pandas (0.20.3)

## Usage

The config file should contain the following information:

```
# EXAMPLE CONFIG FILE FOR MULTI-GENOME OLIGO DESIGN

# target coordinates for oligo design

target_bed: regions.bed

# genome to which the target coordinates correspond

source_genome: hg38.noalt.fa

# GTF for assignment of oligos to strand

gtf: gencode.v32.gtf

# all genomes to screen for homologous sequences (probably this will also include the source genome)

all_genomes:
- hg38.noalt.fa
- panTro6.fa

# shorthand for to each of the genomes for file and sequence names (same order as all_genomes)

names:
- human
- chimp

# filter parameters for doubleFiltBlat.py (for each genome, same order as all_genomes)
# the 6 numbers are...
# first pass min %id, first pass min alignment length (for 200mers), first pass max number of alignments, second pass min %id, second pass min alignment length, second pass max number of alig$

filter_parameters:
- 90 190 4 98 195 4
- 90 190 1 98 195 1

# number of scramble control sequences to generate

n_controls: 1000

# universal priming sequence to add to each end of the oligos

universal_priming_seq_fiveprime: AGGACCGGATCAACT

universal_priming_seq_threeprime: CATTGCGTGAACCGA

# restriction sites to filter out of final oligos

re_sites:
- ACCGGT
- CCTGCAGG
```

The filtering parameters can be adjusted as needed; check the log information generated by `doubleFiltBlat.py`. More information:

```
./doubleFiltBlat.py -h
```

Finally, the three accessory Python scripts should be saved in a directory named `scripts`.

When ready, run the Snakefile:

```
snakemake --configfile config.yaml
```

The final filtered oligos will be savd in the working directory as `final_sequences.fa`!
