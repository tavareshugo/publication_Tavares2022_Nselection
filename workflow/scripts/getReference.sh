#!/bin/bash

#### Input arguments (in order) ####

RELEASE="$1" # ENSEMBL genome release number
OUTDIR="$2"  # output directory name


#### prepare directories ####

mkdir -p "$OUTDIR"
cd "$OUTDIR"


#### Download reference genome ####

# Download genome from ENSEMBL
echo "Downloading reference genome..."
wget -O genome.fa.gz ftp://ftp.ensemblgenomes.org/pub/plants/release-${RELEASE}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Decompress
echo "Decompressing reference fasta..."
gunzip genome.fa.gz


#### Indexing ####

# samtools
echo "Indexing with samtools"
samtools faidx genome.fa

# picard
echo "Indexing with picard"
picard CreateSequenceDictionary REFERENCE=genome.fa OUTPUT=genome.dict

# bwa
echo "Indexing with bwa"
bwa index genome.fa


