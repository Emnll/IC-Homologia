# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 00:53:21 2025

@author: Jean
"""

import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
import os 
from Bio.SeqRecord import SeqRecord

# Função para criação dos arquivos de saída em fasta
def dataframe_to_fasta(df, nome_fasta_file, coluna_id_fasta,coluna_sequencia_fasta, coluna_descricao_fasta):
    records =[]
    with open(nome_fasta_file, 'w') as nome_fasta_file:
        for index, row in df.iterrows():
            record = SeqRecord(Seq(row[coluna_sequencia_fasta]), id=row[coluna_id_fasta],description=row[coluna_descricao_fasta])
            records.append(record)
        SeqIO.write(records,nome_fasta_file,"fasta")
        
# Função que extrai todos os termos entre colchetes []
def extrair_organismos_multiplos(definicao):
    return re.findall(r'\[([^\[\]]+)\]', str(definicao))

# Função para verificar se há algum dos esperados
def tem_organismo_esperado(lista_organismos):
    return any(org in lista_organismos for org in organismos_esperados)
'''
def cria_string_panther(df,gene_column):
    gene_list = df[gene_column].dropna().unique().tolist()
    gene_list_str = ",".join(gene_list)
    return gene_list_str
'''      
path_COG_GO = r"COG_GO\Resultado_cog_GO_left_sem_putative.csv" #Resultado_cog_GO.csv
path_KOG_GO = r"KOG_GO\Resultado_kog_GO_left_sem_putative.csv" #Resultado_kog_GO.csv
path_COSTA_GO = r"COSTA_GO\Resultado_costa_GO_left_sem_putative.csv" #Resultado_costa_GO.csv

diretorio_atual = os.getcwd()

caminho_cog = os.path.join(diretorio_atual,path_COG_GO)
caminho_kog = os.path.join(diretorio_atual,path_KOG_GO)
caminho_costa = os.path.join(diretorio_atual,path_COSTA_GO)

pd_cog = pd.read_csv(caminho_cog,
                         sep=',',
                         low_memory=False
                     )
pd_kog = pd.read_csv(caminho_kog,
                         sep=','
                     )
pd_costa = pd.read_csv(caminho_costa,
                         sep=','  
                       )

# COG e KOG possuem duplicatas para os seus "Protein_ID" e "NCBI_REFSEQ_ACC_VERSION", respectivamente, pois diferentes GenBankID geraram
#    o mesmos valores.
# Costa não há repetições para Code_Gene_DEG porém há: 1318 entradas com Seq_Prot = "Notavailablenow." impossibilitando criação do .fasta destas tuplas.
# Criação com as seguinte qtd de proteías por arquivo:
#    cog: 40556
#    kog: 22346
#    costa: 2967


pd_cog_sem_duplicatas = pd_cog.drop_duplicates(subset="Protein_ID", keep='first')
pd_kog_sem_duplicatas = pd_kog.drop_duplicates(subset="NCBI_REFSEQ_ACC_VERSION", keep='first')
pd_costa_sem_notavailable = pd_costa[pd_costa["Seq_Prot"] != "Notavailablenow."]


#ETAPA COG
#Cog não haverá divisão pois serão analisados como Ortólogos Universais

# entradas: 40551
dataframe_to_fasta(
    pd_cog_sem_duplicatas,
    "cog_pos_go.fasta",
    "Protein_ID",
    "Sequencia",
    "UNIPROTKB_PROTEIN_ACC"
    )

#ETAPA KOG
#[Caenorhabditis elegans],[Homo sapiens],[Drosophila melanogaster], [Saccharomyces cerevisiae S288C]

organismos_esperados = [
    "Caenorhabditis elegans",
    "Homo sapiens",
    "Drosophila melanogaster",
    "Saccharomyces cerevisiae S288C"
]

# Nova coluna com todos os organismos extraídos de cada string
pd_kog_sem_duplicatas["Organismos_extraidos"] = pd_kog_sem_duplicatas["NCBI_DEFINITION"].apply(extrair_organismos_multiplos)

# Cria uma coluna booleana
pd_kog_sem_duplicatas["Tem_organismo_esperado"] = pd_kog_sem_duplicatas["Organismos_extraidos"].apply(tem_organismo_esperado)

# Filtrar:
df_esperados = pd_kog_sem_duplicatas[pd_kog_sem_duplicatas["Tem_organismo_esperado"] == True]

#Todos são de Rattus norvegicus e Mus musculus, gerados a partir de antigos gi humanos que passaram para refseq de ratos. 49 entradas
df_nao_esperados = pd_kog_sem_duplicatas[pd_kog_sem_duplicatas["Tem_organismo_esperado"] == False]

df_human = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Homo sapiens" in x)]
df_droso = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Drosophila melanogaster" in x)]
df_yeast = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Saccharomyces cerevisiae S288C" in x)]
df_elegans = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Caenorhabditis elegans" in x)]



# entradas: 9579
dataframe_to_fasta(
    df_human,
    "kog_humano_pos_go.fasta",
    "NCBI_REFSEQ_ACC_VERSION",
    "NCBI_SEQUENCE",
    "UNIPROTKB_PROTEIN_ACC"
    )
# entradas: 2148
dataframe_to_fasta(
    df_droso,
    "kog_droso_pos_go.fasta",
    "NCBI_REFSEQ_ACC_VERSION",
    "NCBI_SEQUENCE",
    "UNIPROTKB_PROTEIN_ACC"
    )
# entradas: 3548
dataframe_to_fasta(
    df_yeast,
    "kog_cerevisiae_pos_go.fasta",
    "NCBI_REFSEQ_ACC_VERSION",
    "NCBI_SEQUENCE",
    "UNIPROTKB_PROTEIN_ACC"
    )
# entradas: 6106
dataframe_to_fasta(
    df_elegans,
    "kog_elegans_pos_go.fasta",
    "NCBI_REFSEQ_ACC_VERSION",
    "NCBI_SEQUENCE",
    "UNIPROTKB_PROTEIN_ACC"
    )


#ETAPA COSTA
organismos_unicos = pd_costa_sem_notavailable["Organism"].unique()
dfs_por_organismo = {
    org: pd_costa_sem_notavailable[pd_costa_sem_notavailable["Organism"] == org].copy()
    for org in organismos_unicos
}

df_costa_human = dfs_por_organismo.get("Homo sapiens")
df_costa_drosophila = dfs_por_organismo.get("Drosophila melanogaster")
df_costa_yeast = dfs_por_organismo.get("Saccharomyces cerevisiae")

# entradas: 3
pd_costa.columns
dataframe_to_fasta(
    df_costa_human,
    "costa_human_pos_go.fasta",
    "Code_Gene_DEG",
    "Seq_Prot",
    "Uniprot"
    )
# entradas: 206
pd_costa.columns
dataframe_to_fasta(
    df_costa_drosophila,
    "costa_droso_pos_go.fasta",
    "Code_Gene_DEG",
    "Seq_Prot",
    "Uniprot"
    )
# entradas: 1059
pd_costa.columns
dataframe_to_fasta(
    df_costa_yeast,
    "costa_cerevisiae_pos_go.fasta",
    "Code_Gene_DEG",
    "Seq_Prot",
    "Uniprot"
    )

'''
nome_arquivo_kog_human = os.path.join(diretorio_atual,"Resultado_kog_human_para_enrichment.csv")
nome_arquivo_kog_droso = os.path.join(diretorio_atual,"Resultado_kog_droso_para_enrichment.csv")
nome_arquivo_kog_yeast = os.path.join(diretorio_atual,"Resultado_kog_yeast_para_enrichment.csv")
nome_arquivo_kog_elegans = os.path.join(diretorio_atual,"Resultado_kog_elegans_para_enrichment.csv")
nome_arquivo_costa_human = os.path.join(diretorio_atual,"Resultado_costa_human_para_enrichment.csv")
nome_arquivo_costa_droso = os.path.join(diretorio_atual,"Resultado_costa_droso_para_enrichment.csv")
nome_arquivo_costa_yeast = os.path.join(diretorio_atual,"Resultado_costa_yeast_para_enrichment.csv")


df_human.to_csv(nome_arquivo_kog_human,index=False)
df_droso.to_csv(nome_arquivo_kog_droso,index=False)
df_yeast.to_csv(nome_arquivo_kog_yeast,index=False)
df_elegans.to_csv(nome_arquivo_kog_elegans,index=False)
df_costa_human.to_csv(nome_arquivo_costa_human,index=False)
df_costa_drosophila.to_csv(nome_arquivo_costa_droso,index=False)
df_costa_yeast.to_csv(nome_arquivo_costa_yeast,index=False)
'''

'''
string_panther_kog_humano = cria_string_panther(df_human,"UNIPROTKB_PROTEIN_ACC")
string_panther_kog_droso = cria_string_panther(df_droso,"UNIPROTKB_PROTEIN_ACC")
string_panther_kog_yeast = cria_string_panther(df_yeast,"UNIPROTKB_PROTEIN_ACC")
string_panther_kog_elegans = cria_string_panther(df_elegans,"UNIPROTKB_PROTEIN_ACC")
string_panther_costa_humano = cria_string_panther(df_costa_human,"Uniprot")
string_panther_costa_droso = cria_string_panther(df_costa_drosophila,"Uniprot")
string_panther_costa_yeast = cria_string_panther(df_costa_yeast,"Uniprot")
'''