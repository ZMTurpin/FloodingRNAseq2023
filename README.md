# FloodingRNAseq2023
Code written by ZMTurpin in June2023. Trims, aligns, filters, and calculates DEGs from fastq input.

  trim.sh
    trims contaminating adapter sequences from demultiplexed fastq files with cutadapt

  Rsem_align.sh
    Create STAR indices and align to transcript models in RSEM. Estimate transcript- and gene- level mRNA abundances (TPM)

  rpm_RNA_GenomeCov.sh
    #This script prepares rpm-normalized per-base coverage files (bw) for each replicate of an RNA-seq experiment AND merges the replicates.
    #inputs are STAR-aligned read pairs
    #ZMT_JULY27_2023

  filter.log
    fix header in RSEM output tables and filter rows with (log2FoldChange)^2 >= 4 AND padj =< 0.05

  DEseq2.log
    perform differential expression analysis in DEseq2

Code written by ZMTurpin in March 2024. 
  joinToHeatmap.log
    from maizeMine DEG ortholog pairs, build up table of shared gene identifiers, lfc values, and padj values

  logicClustser.py
    assign logical group identification to maize_ortholog_degTable.tsv
    for each genotype, specify each of 3 possible DEG outcomes (up, down, notsignificant) and assign each B73 ortholog to one of these categories
    sort df by logicCat, log2FoldChange_B73, log2FoldChange_Ky21, log2FoldChange_Oh7b, and plot a heatmap of the log2foldchange data
