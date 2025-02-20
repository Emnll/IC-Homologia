# -*- coding: utf-8 -*-
"""
Recebe o arquivo GeneCount e Orthogroups e então encontra os determinados
Orthogroups que possuem ortologia a todos os organismos
por fim retorna arquivos com as proteínas ortólogas a todos os organismos da execução
"""
import pandas as pd
import os

dict_final = {}
diretorio_GeneCount = r"GeneCount"
diretorio_Orthogroups = r"Orthogroups"
diretorio_Merge = r"ResultadosMerge"

diretorio_atual = os.getcwd()
diretorio = os.path.join(diretorio_atual, diretorio_GeneCount)
caminho_Orthogroups = os.path.join(diretorio_atual, diretorio_Orthogroups)
caminho_Merge = os.path.join(diretorio_atual, diretorio_Merge)

def mergeOrtho(df_GeneCount_filtrado, nomeEspecies):
    arquivo = nomeEspecies + ".tsv"
    caminho_Orthogroups_merge_final = os.path.join(caminho_Merge, arquivo)
    
    nome_arquivo_Orthogroups = "Orthogroups." + arquivo
    caminho_arquivo_Orthogroups = os.path.join(caminho_Orthogroups, nome_arquivo_Orthogroups)
    
    df_Orthogroups = pd.read_table(caminho_arquivo_Orthogroups)
    
    df_merge_Ortologos = pd.merge(df_GeneCount_filtrado, df_Orthogroups, on = 'Orthogroup', how = 'left')
    
    df_merge_Ortologos.to_csv(caminho_Orthogroups_merge_final, sep="\t",index=False)
    
    
    # Carregar o arquivo TSV em um DataFrame
for arquivo in os.listdir(diretorio):
   if arquivo.endswith('.tsv'):
       caminho_arquivo = os.path.join(diretorio, arquivo)
       df_original_GeneCount = pd.read_table(caminho_arquivo)
       df_GeneCount = df_original_GeneCount.iloc[:, 1:-1]
       quantidade_Especies = df_GeneCount.shape[1]
       grupoOrt = (df_GeneCount != 0).sum(axis = 1)
       df_GeneCount_filtrado = df_original_GeneCount[grupoOrt == quantidade_Especies]
       contagem_grupoOrt = (df_GeneCount_filtrado['protein_mansoni']).sum()
       key_name = arquivo.split('.')[2]
       dict_final[key_name] = contagem_grupoOrt
       mergeOrtho(df_GeneCount_filtrado.iloc[:, [0,-1]], key_name)
            
pd_final = pd.DataFrame([dict_final])

caminho_Orthogroups_merge_final = os.path.join(caminho_Merge, "tabela_final")
pd_final.to_csv(caminho_Orthogroups_merge_final, sep="\t",index=False)

