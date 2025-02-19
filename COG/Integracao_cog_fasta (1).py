import pandas as pd
import numpy as np
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

diretorio_atual = os.getcwd()
diretorio_cogfa = os.path.join(diretorio_atual, r"cog-20.fa")
diretorio_cogcsv = os.path.join(diretorio_atual, r"cog-20.cog.csv")


df_cog =  pd.read_csv(diretorio_cogcsv, sep=",")
df_cog.columns = ["Gene_ID","NCBI_Ass_ID","Protein_ID","Protein_length","COG_footprint_coordinates","COG_footprint_length","COG_ID","reserved","COG_member_class","PSI-BLAST_score","PSI-BLAST_evalue","COG_profile_length","Protein_footprint_coordinates"]

url_cog_fast = SeqIO.parse(open(diretorio_cogfa),'fasta')




ortologos_universais = np.array((["COG0012",
"COG0016",
"COG0018",
"COG0048",
"COG0049",
"COG0052",
"COG0060",
"COG0080",
"COG0081",
"COG0085",
"COG0087",
"COG0091",
"COG0092",
"COG0093",
"COG0094",
"COG0096",
"COG0097",
"COG0098",
"COG0099",
"COG0100",
"COG0102",
"COG0103",
"COG0124",
"COG0143",
"COG0172",
"COG0184",
"COG0186",
"COG0197",
"COG0200",
"COG0201",
"COG0202",
"COG0256",
"COG0495",
"COG0522",
"COG0525",
"COG0533"]))

df_cog_ortologos_universais = df_cog[df_cog['COG_ID'].isin(ortologos_universais)].copy()
df_cog_ortologos_universais['Protein_ID'] = df_cog_ortologos_universais['Protein_ID'].str.replace('.', '_', regex=False)

fasta_data = [
    {"Protein_ID": record.id, "COG_DESCRICAO": record.description, "COG_SEQUENCE": str(record.seq)}
    for record in url_cog_fast
]
data_proteina = pd.DataFrame(fasta_data)

fun_20_tab = pd.read_table(os.path.join(diretorio_atual,r"fun-20.tab.txt"))
fun_20_tab.columns = ["COG_FUNC","Hex","COG_FUNC_CATEG_DESCR"]
fun_20_tab.drop(["Hex"], axis = 1)

cog_20_def_tab = pd.read_table(os.path.join(diretorio_atual,r"cog-20.def.tab.txt"),encoding='ISO-8859-1')
cog_20_def_tab.columns = ["COG_ID","COG_FUNC","COG_NAME","GENE_ASSOC","FUNC_PATHWAY_COG","PUBMED_ID","PDB_ID"]


merge_cog_universal_fasta = pd.merge(df_cog_ortologos_universais, data_proteina, on = 'Protein_ID', how = 'left')


dataframe_resultante = pd.merge(merge_cog_universal_fasta,cog_20_def_tab, on='COG_ID',how = 'left' )
dataframe_final = pd.merge(dataframe_resultante, fun_20_tab, on='COG_FUNC',how = 'left' )
dataframe_final.to_csv(r"ortologos_universais.csv", sep="\t", index=False, header=False)

output_fasta_path = os.path.join(diretorio_atual, "ortologos_universais.fasta")


fasta_records = [
    SeqRecord(
        Seq(row['COG_SEQUENCE']),
        id=row['Protein_ID'],     
        description=row['COG_DESCRICAO']
    )
    for _, row in dataframe_final.iterrows()
]

# Write to a FASTA file
with open(output_fasta_path, "w") as fasta_file:
    SeqIO.write(fasta_records, fasta_file, "fasta")


