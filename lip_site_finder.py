# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 10:06:48 2023

@author: lawashburn
"""

import csv
import pandas as pd
import re
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
from bioservices import UniProt

db_path = r"C:\Users\lawashburn\Documents\collaborations\Haiyan\phospho_glyco_20231009\uniprot-download_true_format_fasta_query__28Homo_20sapiens_29_20AND_-2022.11.22-18.39.02.76.fasta" #database fasta path
data_path = r"C:\Users\lawashburn\Documents\collaborations\Haiyan\phospho_glyco_20231009\example_input.csv"
output_path = r"C:\Users\lawashburn\Documents\collaborations\Haiyan\phospho_glyco_20231009\output"
PTM_site_variance = 0

fasta_to_df = []
title_to_df = []
accession_to_df = []
with open(db_path) as fasta_file:  # Will close handle cleanly
    for title, sequence in SimpleFastaParser(fasta_file):
        fasta_to_df.append(sequence)
        title_to_df.append(title)


fasta_df = pd.DataFrame()
fasta_df['Record ID'] = title_to_df
fasta_df['Protein Sequence'] = fasta_to_df
fasta_df['Record ID']= fasta_df['Record ID'].str.split('|').str[1].str.strip('|')
#%%
data = pd.read_csv(data_path)
merged = data.merge(fasta_df, left_on=('Leading razor protein_P'),right_on='Record ID', how='inner')
queries = len(merged)
protein_seqs = merged['Protein Sequence'].values.tolist()
peptide_seqs = merged['Sequence S'].values.tolist()
#%%

pep_storage = []
prot_storage = []
start_index = []
end_index = []
peptide_length = []
protein_cleavage_before = []
peptide_first_cleavage = []
peptide_last_cleavage = []
protein_cleavage_after = []
first_cleavage_type = []
second_cleavage_type = []

for a in range(0,queries):
    protein_query = protein_seqs[a]
    peptide_query = peptide_seqs[a]
    
    if peptide_query in protein_query:
        prot_index = protein_query.index(peptide_query)
        start_index.append(prot_index+1)
        pep_len = len(peptide_query)
        peptide_length.append(pep_len)
        end_index.append(prot_index+pep_len)
        
        if prot_index == 0:
            prot_res_before = 'None'
            protein_cleavage_before.append(prot_res_before)
            first_cleavage_type.append('Non-tryptic')
        else:
            prot_res_before = protein_query[prot_index-1]
            protein_cleavage_before.append(prot_res_before)
            
            if prot_res_before == 'K':
                first_cleavage_type.append('Tryptic')
            elif prot_res_before == 'R':
                first_cleavage_type.append('Tryptic')
            else:
                first_cleavage_type.append('Non-tryptic')
            
        first_pep_res = peptide_query[0]
        last_pep_res = peptide_query[pep_len-1]
        
        if last_pep_res == 'K':
            second_cleavage_type.append('Tryptic')
        elif last_pep_res == 'R':
            second_cleavage_type.append('Tryptic')
        else:
            second_cleavage_type.append('Non-tryptic')
        
        if (prot_index+pep_len)==len(protein_query):
            prot_res_after = 'None'
            protein_cleavage_after.append(prot_res_after)
        else:
            prot_res_after = protein_query[prot_index+pep_len]
            protein_cleavage_after.append(prot_res_after)

        peptide_first_cleavage.append(first_pep_res)
        peptide_last_cleavage.append(last_pep_res)
        pep_storage.append(peptide_query)
        prot_storage.append(protein_query)

index_report = pd.DataFrame()
index_report['Peptide sequence'] = pep_storage
index_report['Protein sequence'] = prot_storage
index_report['Peptide start index'] = start_index
index_report['Peptide end index'] = end_index
index_report['Peptide length'] = peptide_length
index_report['Protein cleavage before'] = protein_cleavage_before
index_report['Peptide first residue'] = peptide_first_cleavage
index_report['Peptide last residue'] = peptide_last_cleavage
index_report['Protein cleavage after'] = protein_cleavage_after
index_report['First cleavage site'] = first_cleavage_type
index_report['Second cleavage site'] = second_cleavage_type

conditions = [
    index_report['First cleavage site'].eq('Tryptic') & index_report['Second cleavage site'].eq('Tryptic'),
    index_report['First cleavage site'].eq('Tryptic') & index_report['Protein cleavage after'].eq('None'),
    index_report['Second cleavage site'].eq('Tryptic') & index_report['Protein cleavage before'].eq('None')
]

choices = ['Fully tryptic','Fully tryptic','Fully tryptic']

index_report['Peptide status'] = np.select(conditions, choices, default='Half tryptic')

merged_final = merged.merge(index_report,left_on=['Sequence S','Protein Sequence'],right_on=['Peptide sequence','Protein sequence'])

half_tryp_filter = merged_final.loc[merged_final['Peptide status'] == 'Half tryptic']

begin_nontryp_filter = half_tryp_filter.loc[half_tryp_filter['First cleavage site'] == 'Non-tryptic']
begin_nontryp_filter['LiP site'] = begin_nontryp_filter['Peptide start index']

end_nontryp_filter = half_tryp_filter.loc[half_tryp_filter['Second cleavage site'] == 'Non-tryptic']
end_nontryp_filter['LiP site'] = end_nontryp_filter['Peptide end index'] +1

lip_site_table = pd.concat([begin_nontryp_filter,end_nontryp_filter])

applicable_mod_storage = []

for index, row in lip_site_table.iterrows():
    prot_query = row['Leading razor protein_P']
    index_query = row['LiP site']

    u = UniProt()
    res = u.get_df(prot_query.split())
    ptm = res['Modified residue'].iloc[0]
    ptm_split = ptm.split('MOD_RES ')
    ptm_split = ptm_split[1:]
    ptm_df = pd.DataFrame(
             {'Res #': ptm_split})
    ptm_df[['Res #','type','evidence']] = ptm_df['Res #'].str.split('; /',expand=True)


    ptm_df['type'] = ptm_df['type'].replace({'note=': '','"':''}, regex=True)

    res_begin = index_query - PTM_site_variance
    res_end = index_query + PTM_site_variance

    ptm_df["Res #"] = pd.to_numeric(ptm_df["Res #"])
    ptm_df_filtered = ptm_df.loc[ptm_df['Res #'] >= res_begin]
    ptm_df_filtered = ptm_df_filtered.loc[ptm_df['Res #'] <= res_end]
    ptm_df_filtered['IDed sites'] = ptm_df_filtered['type'] + '(' + ptm_df_filtered['Res #'].astype(str) + ')'
    
    applicable_mods = ptm_df_filtered['IDed sites'].values.tolist()
    
    if len(applicable_mods)>0:
        applicable_mod_storage.append(applicable_mods)
    else:
        applicable_mod_storage.append('No modifications found')

lip_site_table['Modifications'] = applicable_mod_storage

lip_site_table = lip_site_table.drop(columns=['Protein sequence','Peptide sequence','Peptide start index','Peptide end index','Peptide length',
                             'Protein cleavage before','Peptide first residue','Peptide last residue','Protein cleavage after','Protein Sequence'])

file_out_path = output_path + '\\lip_site_PTM_finder_results.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        lip_site_table.to_csv(filec,index=False)

        