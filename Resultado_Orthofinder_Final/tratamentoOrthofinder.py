# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 20:54:16 2024

@author: Emanu
"""
#%%
import pandas as pd
import os

#%%
""" Dicionarios """
dict_final = {}
dict_final_hsa = {}
visualizer = pd.DataFrame()

#%%
""" Trocar pasta """
diretorio_pasta = r"OU"

#%%
""" Diretorios das pastas """
diretorio_GeneCount = r"GeneCount"
diretorio_Orthogroups = r"Orthogroups"
diretorio_Merge = r"ResultadosMerge"
diretorio_atual = os.path.join(os.getcwd(), diretorio_pasta)

#%%
""" Caminhos dos diretorios """
diretorio = os.path.join(diretorio_atual, diretorio_GeneCount)
caminho_Orthogroups = os.path.join(diretorio_atual, diretorio_Orthogroups)
caminho_Merge = os.path.join(diretorio_atual, diretorio_Merge)

#%%
def mergeOrtho(df_GeneCount_filtrado, nomeEspecies):
    arquivo = nomeEspecies + ".tsv"
    caminho_Orthogroups_merge_final = os.path.join(caminho_Merge, arquivo)
    
    nome_arquivo_Orthogroups = "Orthogroups." + arquivo
    caminho_arquivo_Orthogroups = os.path.join(caminho_Orthogroups, nome_arquivo_Orthogroups)
    
    df_Orthogroups = pd.read_table(caminho_arquivo_Orthogroups)
    
    df_merge_Ortologos = pd.merge(df_GeneCount_filtrado, df_Orthogroups, on = 'Orthogroup', how = 'left')
    
    df_merge_Ortologos.to_csv(caminho_Orthogroups_merge_final, sep="\t",index=False)

#%%
def flatten_extend(matrix):
    flat_list = []
    for sublist in matrix:
        if isinstance(sublist, list):
            flat_list.extend(sublist)
        else:
            flat_list.append(sublist)
    return flat_list
    
#%%    
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

#%%            
pd_final = pd.DataFrame([dict_final])

#%%
caminho_Orthogroups_merge_final = os.path.join(caminho_Merge, "tabela_final")
pd_final.to_csv(caminho_Orthogroups_merge_final, sep="\t",index=False)

#%%
for arquivo_merge in os.listdir(caminho_Merge):
    nome = arquivo_merge.strip('.tsv')
    if arquivo_merge.endswith('.tsv'):
        dict_final_hsa[nome] = []
        caminho_arquivo = os.path.join(caminho_Merge, arquivo_merge)
        visualizer = pd.read_table(caminho_arquivo)
        teste = []
        for k, v in visualizer.iterrows():
            if "protein_mansoni" in v:
                teste.append(v['protein_mansoni'].split(', '))
            teste = flatten_extend(teste)
            dict_final_hsa[nome] = teste

#%%
""" Tirar hsa """
#Obs: Não funciona com os OU
hsa_df = pd.DataFrame(dict_final_hsa['manxhsa'], columns = ['Proteins'])
del dict_final_hsa['manxhsa']

#%%
""" """

for key, values in dict_final_hsa.items():
    df = pd.DataFrame(dict_final_hsa[key], columns = ['Proteins'])
    mergido = df.merge(hsa_df, on = ['Proteins'], how = 'left', indicator = True)
    menos_hsa = mergido[mergido['_merge'] == 'left_only'].drop(columns = ['_merge'])

    filename = f"{diretorio_atual}/Resultado_menos_hsa_essencial/{key}_resultado_final.txt"

    menos_hsa.to_csv(filename, sep="\t", index=False, header=False)
    print(f"Created {filename}")
# %%
