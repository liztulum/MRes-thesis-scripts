# -*- coding: utf-8 -*-
"""
Spyder Editor
 
This is a temporary script file.
"""
 
import json
import pandas as pd
 
# add your file path to the json file
gene_file = open('C:/Users/liz.tulum/OneDrive - Unilever/MRes course/190620_Human_Whole_Transcriptome_2.0_Manifest_probe_to_gene (3).json')
 
# load in json file of gene names as data frame
gene_names = json.load(gene_file)
 
# the concentrations of the file your looking add
niac_conc= ["0.512", "2.56", "12.8", "64", "320", "1600", "8000"]
niac_conc2= ["3.84", "19.2", "96", "480", "2400", "12000", "60000"]
doxo_conc= ["6.4e.05", "0.00032", "0.0016", "0.008", "0.04", "0.2", "1"]
 
# for each concentration it reads in the file as a data frame and replaces the probe name with the gene name
# HepaRG
for conc in niac_conc:
    df = pd.read_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepaRG_Niacinamide_DESeq_CONCENTRATION_{conc}_vs_0.txt', sep='\t')
    print(df)
    # change 'Unnamed:0' to 'gene' or whatever the first column of the file is caaled
    df['gene'].replace(gene_names,inplace = True)
   
    # export df with the gene names as a tsv
    df.to_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepaRG_Niacinamide_DESeq_CONCENTRATION_{conc}_vs_0.csv', index=False)
 
for conc in doxo_conc:
    df = pd.read_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepaRG_Doxorubicin_HCl_DESeq_CONCENTRATION_{conc}_vs_0.txt', sep='\t')
    # change 'Unnamed:0' to 'gene' or whatever the first column of the file is caaled
    df['gene'].replace(gene_names,inplace = True)
    # export df with the gene names as a tsv
    df.to_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepaRG_Doxorubicin_HCl_DESeq_CONCENTRATION_{conc}_vs_0.csv', index=False)
 
#HepG2
for conc in niac_conc2:
    df = pd.read_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepG2_Niacinamide_DESeq_CONCENTRATION_{conc}_vs_0.txt', sep='\t')
    # change 'Unnamed:0' to 'gene' or whatever the first column of the file is caaled
    df['gene'].replace(gene_names,inplace = True)
   
    # export df with the gene names as a tsv
    df.to_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepG2_Niacinamide_DESeq_CONCENTRATION_{conc}_vs_0.csv', index=False)
 
for conc in doxo_conc:
    df = pd.read_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepG2_Doxorubicin_HCl_DESeq_CONCENTRATION_{conc}_vs_0.txt', sep='\t')
    # change 'Unnamed:0' to 'gene' or whatever the first column of the file is caaled
    df['gene'].replace(gene_names,inplace = True)
    # export df with the gene names as a tsv
    df.to_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/HepG2_Doxorubicin_HCl_DESeq_CONCENTRATION_{conc}_vs_0.csv', index=False)
 
#MCF7
for conc in niac_conc2:
    df = pd.read_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/MCF-7_Niacinamide_DESeq_CONCENTRATION_{conc}_vs_0.txt', sep='\t')
    # change 'Unnamed:0' to 'gene' or whatever the first column of the file is caaled
    df['gene'].replace(gene_names,inplace = True)
   
    # export df with the gene names as a tsv
    df.to_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/MCF-7_Niacinamide_DESeq_CONCENTRATION_{conc}_vs_0.csv', index=False)
 
for conc in doxo_conc:
    df = pd.read_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/MCF-7_Doxorubicin_HCl_DESeq_CONCENTRATION_{conc}_vs_0.txt', sep='\t')
    # change 'Unnamed:0' to 'gene' or whatever the first column of the file is caaled
    df['gene'].replace(gene_names,inplace = True)
    # export df with the gene names as a tsv
    df.to_csv(f'C:/Users/liz.tulum/OneDrive - Unilever/MRes course/MCF-7_Doxorubicin_HCl_DESeq_CONCENTRATION_{conc}_vs_0.csv', index=False)
