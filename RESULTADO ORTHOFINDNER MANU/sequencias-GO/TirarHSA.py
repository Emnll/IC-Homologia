# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 09:58:29 2025

@author: Emanu
"""

import pandas as pd
import os

dict_final = {}
visualizer = pd.DataFrame()

diretorio_pasta = r"cog_pos_go"
diretorio_GeneCount = r"GeneCount"
diretorio_Orthogroups = r"Orthogroups"
diretorio_Merge = r"ResultadosMerge"
diretorio_atual = os.path.join(os.getcwd(), diretorio_pasta)

diretorio = os.path.join(diretorio_atual, diretorio_GeneCount)
caminho_Orthogroups = os.path.join(diretorio_atual, diretorio_Orthogroups)
caminho_Merge = os.path.join(diretorio_atual, diretorio_Merge)

def flatten_extend(matrix):
    """Flatten a list of lists into a single list."""
    flat_list = []
    for sublist in matrix:
        if isinstance(sublist, list):  # Ensure we only flatten lists
            flat_list.extend(sublist)
        else:  # If it's not a list, append as is
            flat_list.append(sublist)
    return flat_list
 
for arquivo_merge in os.listdir(caminho_Merge):
    nome = arquivo_merge.strip('.tsv')
    if arquivo_merge.endswith('.tsv'):
        dict_final[nome] = []
        caminho_arquivo = os.path.join(caminho_Merge, arquivo_merge)
        visualizer = pd.read_table(caminho_arquivo)
        teste = []
        for k, v in visualizer.iterrows():
            if "protein_mansoni" in v:
                teste.append(v['protein_mansoni'].split(', '))
            teste = flatten_extend(teste)
            dict_final[nome] = teste
            
hsa_df = pd.DataFrame(dict_final['manxhsa'], columns = ['Proteins'])
hsa_df2 = pd.Data()
del dict_final['manxhsa']

for key, values in dict_final.items():
    df = pd.DataFrame(dict_final[key], columns = ['Proteins'])
    mergido = df.merge(hsa_df, on = ['Proteins'], how = 'left', indicator = True)
    menos_hsa = mergido[mergido['_merge'] == 'left_only'].drop(columns = ['_merge'])
    filename = f"Resultado_menos_hsa/{key}_resultado_final.txt"
    menos_hsa.to_csv(filename, sep="\t", index=False, header=False)
    print(f"Created {filename}")