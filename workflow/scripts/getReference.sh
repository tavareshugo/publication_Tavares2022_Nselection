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


#### Centromere annotation ####

# Centromere region is annotated in TAIR9 release
# see https://www.biostars.org/p/18782/
wget -O temp_centromeres.txt ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_gff3/Assembly_GFF/TAIR9_GFF3_assemblies.gff

# extract centromere annotations
printf "chrom\tstart\tend\n" > centromeres.tsv
grep "CEN" temp_centromeres.txt | cut -f 1,4,5 >> centromeres.tsv
rm temp_centromeres.txt
