# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 15:07:51 2025

@author: Jean
"""
#%%

import os
import re
import pandas as pd
import json
#import networkx as nx
#import matplotlib.pyplot as plt
#from networkx.drawing.nx_agraph import graphviz_layout
import logging
from itertools import product

#%%
def extrair_ids_unicos(df, colunas):
    """Extrai IDs únicos de uma ou mais colunas de um dataframe."""
    ids = set()
    # Só pega colunas que realmente existem
    colunas = [c for c in colunas if c in df.columns]
    for coluna in colunas:
        # dropna + astype(str) para garantir strings
        for valores in df[coluna].dropna().astype(str):
            # separa por vírgula; troque por regex se tiver separadores mistos
            for v in valores.split(","):
                v = v.strip()
                if v:
                    ids.add(v)
    return ids

#%%
def contem_id(x, ids_unicos):
    """Checa se x (string com ids separados por vírgula) contém algum id de ids_unicos."""
    if pd.isna(x):
        return False
    # separa e normaliza
    valores = set(s.strip() for s in str(x).split(",") if s.strip())
    return not ids_unicos.isdisjoint(valores)

#%%
#Encontrar interseção entre proteínas de uma mesma etapa porém com retirada de humanos distintos
def contem_proteina(x,proteinas_intersecao_set):
    valores = str(x).split(",")  # ou outro separador
    return any(proteina in valores for proteina in proteinas_intersecao_set)


#%%
pasta_base = os.getcwd()

""" Trocar pasta """
diretorio_pasta = r"Jessica"
diretorio_pasta2 = r"Kog"
diretorio_pasta3 = r"OU"

""" Diretorios das pastas """
diretorio_Resultado_menos_hsa_essencial = r"Resultado_menos_hsa_essencial"
diretorio_Resultado_menos_hsa_todo = r"Resultado_menos_hsa_todo"


diretorio_atual = os.path.join(os.getcwd(), diretorio_pasta)
diretorio_atual2 = os.path.join(os.getcwd(), diretorio_pasta2)
diretorio_atual3 = os.path.join(os.getcwd(), diretorio_pasta3)

""" Caminhos dos diretorios """
diretorio_jessica1 = os.path.join(diretorio_atual, diretorio_Resultado_menos_hsa_essencial)
diretorio_jessica2 = os.path.join(diretorio_atual, diretorio_Resultado_menos_hsa_todo)

diretorio_kog1 = os.path.join(diretorio_atual2, diretorio_Resultado_menos_hsa_essencial)
diretorio_kog2 = os.path.join(diretorio_atual2, diretorio_Resultado_menos_hsa_todo)

diretorio_ou1 = os.path.join(diretorio_atual3, diretorio_Resultado_menos_hsa_essencial)
diretorio_ou2 = os.path.join(diretorio_atual3, diretorio_Resultado_menos_hsa_todo)


#%%
#df = pd.read_csv( os.path.join(caminho,"Todo_resultado_final.txt"),sep="\t", names=['Proteina'])
#print(df.nunique())

df1 = {}
df2 = {}
df3 = {}
df4 = {}
df5 = {}
df6 = {}

for arquivo in os.listdir(diretorio_ou1):
   p = pd.read_csv( os.path.join(diretorio_ou1,arquivo),sep="\t", names=['Proteina'])
   df5[arquivo] = p
   print(f"Quantidade {diretorio_ou1} $$$$ {arquivo} PRIMEIRO: ",p.nunique())

for arquivo in os.listdir(diretorio_ou2):
   p = pd.read_csv( os.path.join(diretorio_ou2,arquivo),sep="\t", names=['Proteina'])
   df6[arquivo] = p
   print(f"Quantidade {diretorio_ou2} $$$$ {arquivo} SEGUNDO: ",p.nunique())
   
for arquivo in os.listdir(diretorio_jessica1):
   p = pd.read_csv( os.path.join(diretorio_jessica1,arquivo),sep="\t", names=['Proteina'])
   df1[arquivo] = p
   print(f"Quantidade {diretorio_jessica1} $$$$ {arquivo} PRIMEIRO: ",p.nunique())

for arquivo in os.listdir(diretorio_jessica2):
    p = pd.read_csv( os.path.join(diretorio_jessica2,arquivo),sep="\t", names=['Proteina'])
    df2[arquivo] = p
    print(f"Quantidade {diretorio_jessica2} $$$$ {arquivo} SEGUNDO: ",p.nunique())
   
for arquivo in os.listdir(diretorio_kog1):
    p = pd.read_csv( os.path.join(diretorio_kog1,arquivo),sep="\t", names=['Proteina'])
    df3[arquivo] = p
    print(f"Quantidade {diretorio_kog1} $$$$ {arquivo} TERCEIRO: ",p.nunique())
   
for arquivo in os.listdir(diretorio_kog2):
    p = pd.read_csv( os.path.join(diretorio_kog2,arquivo),sep="\t", names=['Proteina'])
    df4[arquivo] = p
    print(f"Quantidade {diretorio_kog2} $$$$ {arquivo} QUARTO: ",p.nunique())
   
   
#df's jessica
#DF1 -> costa - hsa_costa // dme => 359; sce => 877; todos => 250
#DF2 -> costa - hsa_kog // dme => 140; sce => 185; todos => 97
#%%
#path pasta merge
pasta_merge_jessica = os.path.join(diretorio_atual,"ResultadosMerge")


dataframe_manxdme_jessica = pd.read_csv(os.path.join(pasta_merge_jessica,"manxdme.tsv"),sep="\t")
dataframe_manxsce_jessica = pd.read_csv(os.path.join(pasta_merge_jessica,"manxsce.tsv"),sep="\t")
dataframe_manxtodos_jessica = pd.read_csv(os.path.join(pasta_merge_jessica,"Todos.tsv"),sep="\t")

for arquivo in os.listdir(pasta_merge_jessica):
    if "manxdme.tsv":
        
        proteinas_df1 = set(df1["manxdme_resultado_final.txt"]["Proteina"])
        proteinas_df2 = set(df2["manxdme_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_dme = proteinas_df1.intersection(proteinas_df2)
        
        df_filtrado_costa_dme = dataframe_manxdme_jessica[dataframe_manxdme_jessica["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_dme))]
        
    if "manxsce.tsv":
        
        proteinas_df1 = set(df1["manxsce_resultado_final.txt"]["Proteina"])
        proteinas_df2 = set(df2["manxsce_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_sce = proteinas_df1.intersection(proteinas_df2)
        
        df_filtrado_costa_sce = dataframe_manxsce_jessica[dataframe_manxsce_jessica["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_sce))]
        
    if "Todos.tsv":
        proteinas_df1 = set(df1["Todo_resultado_final.txt"]["Proteina"])
        proteinas_df2 = set(df2["Todo_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_total = proteinas_df1.intersection(proteinas_df2)
        
        df_filtrado_costa_total = dataframe_manxtodos_jessica[dataframe_manxtodos_jessica["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_total))]
        


#df's kog
#DF3 -> kog - hsa_costa // cel => 2869; dme => 1889; sce => 2274; todos => 860
#DF4 -> kog - hsa_kog // cel => 622; dme => 385; sce => 492; todos => 207
pasta_merge_kog = os.path.join(diretorio_atual2,"ResultadosMerge")


dataframe_manxcel_kog = pd.read_csv(os.path.join(pasta_merge_kog,"manxcel.tsv"),sep="\t")
dataframe_manxdme_kog = pd.read_csv(os.path.join(pasta_merge_kog,"manxdme.tsv"),sep="\t")
dataframe_manxsce_kog = pd.read_csv(os.path.join(pasta_merge_kog,"manxsce.tsv"),sep="\t")
dataframe_manxtodos_kog = pd.read_csv(os.path.join(pasta_merge_kog,"Todos.tsv"),sep="\t")

for arquivo in os.listdir(pasta_merge_kog):
    if "manxcel.tsv":
        proteinas_df3 = set(df3["manxcel_resultado_final.txt"]["Proteina"])
        proteinas_df4 = set(df4["manxcel_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_cel_kog = proteinas_df3.intersection(proteinas_df4)
        
        df_filtrado_kog_cel = dataframe_manxcel_kog[dataframe_manxcel_kog["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_cel_kog))]
    if "manxdme.tsv":
        
        proteinas_df3 = set(df3["manxdme_resultado_final.txt"]["Proteina"])
        proteinas_df4 = set(df4["manxdme_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_dme_kog = proteinas_df3.intersection(proteinas_df4)
        
        df_filtrado_kog_dme = dataframe_manxdme_kog[dataframe_manxdme_kog["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_dme_kog))]
        
    if "manxsce.tsv":
        
        proteinas_df3 = set(df3["manxsce_resultado_final.txt"]["Proteina"])
        proteinas_df4 = set(df4["manxsce_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_sce_kog = proteinas_df3.intersection(proteinas_df4)
        
        df_filtrado_kog_sce = dataframe_manxsce_kog[dataframe_manxsce_kog["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_sce_kog))]
        
    if "Todos.tsv":
        proteinas_df3 = set(df3["Todo_resultado_final.txt"]["Proteina"])
        proteinas_df4 = set(df4["Todo_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_total_kog = proteinas_df3.intersection(proteinas_df4)
        
        df_filtrado_kog_total = dataframe_manxtodos_kog[dataframe_manxtodos_kog["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_total_kog))]
        
#df's OU
#DF5 -> OU - hsa_costa // variedade de espécies primariamente de procariotos
#DF6 -> OU - hsa_kog // variedade de espécies primariamente de procariotos
pasta_merge_ou = os.path.join(diretorio_atual3,"ResultadosMerge")


dataframe_OU = pd.read_csv(os.path.join(pasta_merge_ou,"OU.tsv"),sep="\t")


for arquivo in os.listdir(pasta_merge_ou):
    if "manxcel.tsv":
        proteinas_df5 = set(df5["OU_resultado_final.txt"]["Proteina"])
        proteinas_df6 = set(df6["OU_resultado_final.txt"]["Proteina"])
        
        proteinas_intersecao_OU = proteinas_df5.intersection(proteinas_df6)
        
        df_filtrado_OU = dataframe_OU[dataframe_OU["protein_mansoni"].apply(lambda x: contem_proteina(x,proteinas_intersecao_OU))]


#CONTAGEM qtd de proteínas distintas na interseção de hsa_costa e hsa_kog
dataframes = {
    "costa_dme":df_filtrado_costa_dme,
    "costa_sce":df_filtrado_costa_sce,
    "costa_total":df_filtrado_costa_total,
    "kog_cel":df_filtrado_kog_cel,
    "kog_dme":df_filtrado_kog_dme,
    "kog_sce":df_filtrado_kog_sce,
    "kog_total":df_filtrado_kog_total,
    "OU":df_filtrado_OU
}
dataframes_columns = {
    "costa_dme":["costa_droso_pos_go"],
    "costa_sce":["costa_cerevisiae_pos_go"],
    "costa_total":["costa_cerevisiae_pos_go","costa_droso_pos_go"],
    "kog_cel":["kog_elegans_pos_go"],
    "kog_dme":["kog_droso_pos_go"],
    "kog_sce":["kog_cerevisiae_pos_go"],
    "kog_total":["kog_cerevisiae_pos_go","kog_droso_pos_go","kog_elegans_pos_go"],
    "OU":["cog_pos_go"]
}


resultados_qtd = {}

for nome, df in dataframes.items():
    todos_ids = ",".join(df["protein_mansoni"].dropna().astype(str))
    ids_unicos = set(todos_ids.split(","))
    resultados_qtd[nome] = len(ids_unicos)

print(resultados_qtd)



lista_ids_df_filtrado_costa_dme = (
    df_filtrado_costa_dme["costa_droso_pos_go"]
    .str.split(", ")        
    .explode()              
    .dropna()               
    .tolist()               
)
lista_ids_df_filtrado_costa_sce = (
    df_filtrado_costa_sce["costa_cerevisiae_pos_go"]
    .str.split(", ")        
    .explode()              
    .dropna()               
    .tolist()               
)

lista_ids_df_filtrado_kog_cel = (
    df_filtrado_kog_cel["kog_elegans_pos_go"]
    .str.split(", ")        
    .explode()              
    .dropna()               
    .tolist()               
)

lista_ids_df_filtrado_kog_dme = (
    df_filtrado_kog_dme["kog_droso_pos_go"]
    .str.split(", ")        
    .explode()              
    .dropna()               
    .tolist()               
)

lista_ids_df_filtrado_kog_sce = (
    df_filtrado_kog_sce["kog_cerevisiae_pos_go"]
    .str.split(", ")        
    .explode()              
    .dropna()               
    .tolist()               
)

lista_ids_df_filtrado_cog = (
    df_filtrado_OU["cog_pos_go"]
    .str.split(", ")        
    .explode()              
    .dropna()               
    .tolist()               
)


#precisaria ver pelos dois 
lista_ids_df_filtrado_costa_total = df_filtrado_costa_total["costa_cerevisiae_pos_go"].str.split(", ").explode().tolist() +\
                                    df_filtrado_costa_total["costa_droso_pos_go"].str.split(", ").explode().tolist()

lista_ids_df_filtrado_kog_total = df_filtrado_kog_total["kog_cerevisiae_pos_go"].str.split(", ").explode().tolist() +\
                                    df_filtrado_kog_total["kog_droso_pos_go"].str.split(", ").explode().tolist() +\
                                    df_filtrado_kog_total["kog_elegans_pos_go"].str.split(", ").explode().tolist()
    
    
#
#
#
#                   ARQUIVOS GERADOS PELO CODIGO QUE ESTA EM 'retirada_putative.py'
#
#


#Separar id's proteínas dos organismos modelos para busca no dataset que contém os uniprots

path_origin_kog = r"Resultado_kog_GO_left_sem_putative.csv"
path_origin_cog = r"Resultado_cog_GO_left_sem_putative.csv"
path_origin_costa = r"Resultado_costa_GO_left_sem_putative.csv"         

df_origin_sem_putative_kog   = pd.read_csv(path_origin_kog) 
df_origin_sem_putative_cog   = pd.read_csv(path_origin_cog)     
df_origin_sem_putative_costa = pd.read_csv(path_origin_costa) 

# --- datasets de referência (coluna alvo de cada um) ---
refs = {
    "kog":   (df_origin_sem_putative_kog,   "NCBI_REFSEQ_ACC_VERSION"),
    "costa": (df_origin_sem_putative_costa, "Code_Gene_DEG"),
    "OU":    (df_origin_sem_putative_cog,   "Protein_ID"),
}

# --- onde vamos guardar tudo, SEM sobrescrever ---
resultados = { "kog": {}, "costa": {}, "OU": {} }

# Se quiser limitar quais nomes entram em cada bucket, declare aqui:
grupos = {
    "kog":   {"kog_cel","kog_dme","kog_sce","kog_total"},
    "costa": {"costa_dme","costa_sce","costa_total"},
    "OU":    {"OU"},
}


# --- loop único, organizado por bucket ---
for bucket, (df_ref, col_ref) in refs.items():
    nomes_para_bucket = grupos.get(bucket, set(dataframes.keys()))
    for nome in nomes_para_bucket:
        df = dataframes[nome]
        colunas = dataframes_columns.get(nome, [])
        ids_unicos = extrair_ids_unicos(df, colunas)

        mask = df_ref[col_ref].apply(lambda x: contem_id(x, ids_unicos))
        df_filtrado = df_ref.loc[mask].copy()

        resultados[bucket][nome] = df_filtrado



# Dicionário para guardar os IDs concatenados
ids_concatenados = {}

for grupo, dfs in resultados.items():
    for nome, df in dfs.items():
        # Descobre qual coluna usar (Uniprot ou UNIPROTKB_PROTEIN_ACC)
        if "Uniprot" in df.columns:
            coluna = "Uniprot"
        elif "UNIPROTKB_PROTEIN_ACC" in df.columns:
            coluna = "UNIPROTKB_PROTEIN_ACC"
        else:
            continue 

        ids_str = ",".join(df[coluna].dropna().astype(str))
        ids_concatenados[f"{grupo}_{nome}"] = ids_str

# Agora você pode acessar assim:
print(ids_concatenados["costa_costa_dme"])
print(ids_concatenados["kog_kog_total"])
print(ids_concatenados["OU_OU"])

#BUSCA DOS ARQUIVOS JSON NO PATNTHER PELAS STRINGS FEITAS EM ids_concatenados(MENOS AS STRINGS COM ORGANISMOS MÚLTIPLOS: costa_total, kog_total, OU).
# Disponível em: https://pantherdb.org/
 
#Utilização da função de Statistical overrepresentation test: 
    #Complete GO annotation datasets - Complete GO annotations including both manually curated and electronic annotations. 
    #Electronic annotations are generated by computer algorithm based on sequence similarity, and are usually not reviewed by curators. 
    #The datasets include all three GO aspects:
        #Biological process
        #Cellular Component
        #Molecular Function
        
    #PANTHER GO-Slim - GO annotations from the phylogenetic curation effort captured in about 3000 GO Slim terms. 
    #The annotation data sets use the data from the Gene Ontology phylogenetic annotation effort. 
    #These are manually curated annotations to the ancestral nodes on PANTHER family trees based on the experimental 
    #   annotations on their leaf descendants (extant genes). The annotation can be inferred to other leaf sequences 
    #   (that have not been tested experimentally) under the annotated ancestral node.
    #The annotations datasets are in all three GO aspects:
        
        #Biological process
        #Cellular Component
        #Molecular Function
        
#Carregamento dos arquivos json advindos da busca pela aplicação web.
#pré processamento dos arquivos.
#Arquivos de GO completo, foram nomeados: <#>_go_complete.<json,txt> e foram alterados os campos list_name para conter os ids das proteínas; 
#   Arquivos sem go_complete foi utilizado go_slim.


panther_results = r"PANTHER_RESULTADOS"

caminho_panther_resultados = os.path.dirname(os.getcwd())

diretorio_arquivos_panther = os.path.join(caminho_panther_resultados, panther_results)  

dicionario_analises_arquivos_panther = {}

for dirpath, dirnames, filenames in os.walk(diretorio_arquivos_panther):
    for file in filenames:
        file_path = os.path.join(dirpath, file)
        
        # Quebra o caminho em partes para juntar após e criar o nome do arquivo para o dicionario
        parts = file_path.split(os.sep)
        ultimos_3 = os.path.join(*parts[-3:])
        nome_limpo = str(ultimos_3).replace("\\", "_").replace(" ", "_")
        
        if file.endswith(".json"):
            with open(file_path, 'r') as f:
                loaded_data = json.load(f)
            dicionario_analises_arquivos_panther[nome_limpo] = loaded_data


print(dicionario_analises_arquivos_panther.keys())




panther_results = r"PANTHER_RESULTADOS"

diretorio_arquivos_panther = os.path.join(os.path.dirname(pasta_base),panther_results)#os.path.join(pasta_base, panther_results)  

dicionario_analises_arquivos_panther = {}

for dirpath, dirnames, filenames in os.walk(diretorio_arquivos_panther):
    for file in filenames:
        file_path = os.path.join(dirpath, file)
        
        # Quebra o caminho em partes para juntar após e criar o nome do arquivo para o dicionario
        parts = file_path.split(os.sep)
        ultimos_3 = os.path.join(*parts[-3:])
        nome_limpo = str(ultimos_3).replace("/", "_").replace(" ", "_")
        print("nomel limpo:",nome_limpo) 
        if file.endswith(".json"):
            with open(file_path, 'r') as f:
                loaded_data = json.load(f)
            dicionario_analises_arquivos_panther[nome_limpo] = loaded_data



# Set up logging for debugging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

dataframe_go_panther_complete_slim = pd.DataFrame(columns=[
    "Origem_organismo",
    "Organismo",
    "Annotation_Type",
    "Statistical_Test_Type",
    "Bonferroni_Correction_value",
    "Data_Version",
    "Panther_GO_id",
    "Term_Label",
    "Term_Level",
    "QTD_proteins_in_list",
    "Expected",
    "Fold_Enrichment",
    "pValue",
    "Plus_Minus",
    "QTD_proteins_in_Reference",
    "Uniprot",
])

rows = []
for key, value in dicionario_analises_arquivos_panther.items():
    try:
        name_temp = key.split("\\")
        origem = name_temp[0].split('_')[0]
        
        # Safely access top-level keys with defaults
        overrepresentation = value.get("overrepresentation", {})
        organismo = overrepresentation.get("upload_lists", {}).get("input_list", {}).get("organism", "Unknown")
        Annotation_Type = overrepresentation.get("annotation_type", "Unknown")
        Bonferroni_Correction_value = overrepresentation.get("bonferroni_correction", "Unknown")
        Statistical_Test_Type = overrepresentation.get("test_type", "Unknown")
        Data_Version = overrepresentation.get("data_version_release_date", "Unknown")
        
        # Get group list, default to empty list if missing
        groups = overrepresentation.get("group", [])
        logging.info(f"Processing {len(groups)} groups for key: {key}")
        
        for grupos in groups:
            # Ensure grupos is a dictionary
            if not isinstance(grupos, dict):
                logging.warning(f"Skipping invalid group: {grupos}")
                continue
                
            # Handle both list and dict cases for grupos["result"]
            result = grupos.get("result")
            if result is None:
                logging.warning(f"No 'result' key in group: {grupos}")
                continue
            result_list = result if isinstance(result, list) else [result]
            
            for result_dict in result_list:
                # Ensure result_dict is a dictionary
                if not isinstance(result_dict, dict):
                    logging.warning(f"Skipping invalid result_dict: {result_dict}")
                    continue
                
                # Get term, handle UNCLASSIFIED or missing term
                term = result_dict.get("term")
                if term == {"label": "UNCLASSIFIED"} or term is None:
                    logging.info(f"Skipping UNCLASSIFIED or missing term in result: {result_dict}")
                    continue
                if isinstance(term, dict):
                    Panther_GO_id = term.get("id", "NA")
                    Term_Label = term.get("label", "NA")
                    Term_Level = term.get("level", "NA")
                else:
                    Panther_GO_id = "unclassified"
                    Term_Label = "unclassified"
                    Term_Level = "NA"
                
                # Safely access input_list and its nested keys
                input_list = result_dict.get("input_list", {})
                mapped_id_list = input_list.get("mapped_id_list", {}).get("mapped_id", [])
                Protein_IDs = ",".join(mapped_id_list) if isinstance(mapped_id_list, list) else str(mapped_id_list)
                
                # Append row with all required fields
                rows.append({
                    "Origem_organismo": origem,
                    "Organismo": organismo,
                    "Annotation_Type": Annotation_Type,
                    "Statistical_Test_Type": Statistical_Test_Type,
                    "Bonferroni_Correction_value": Bonferroni_Correction_value,
                    "Data_Version": Data_Version,
                    "Panther_GO_id": Panther_GO_id,
                    "Term_Label": Term_Label,
                    "Term_Level": Term_Level,
                    "QTD_proteins_in_list": input_list.get("number_in_list", "NA"),
                    "Expected": input_list.get("expected", "NA"),
                    "Fold_Enrichment": input_list.get("fold_enrichment", "NA"),
                    "pValue": input_list.get("pValue", "NA"),
                    "Plus_Minus": input_list.get("plus_minus", "NA"),
                    "QTD_proteins_in_Reference": result_dict.get("number_in_reference", "NA"),
                    "Uniprot": Protein_IDs
                })
                
                # Batch DataFrame creation for memory efficiency
                if len(rows) >= 1000:
                    temp_df = pd.DataFrame(rows)
                    dataframe_go_panther_complete_slim = pd.concat([dataframe_go_panther_complete_slim, temp_df], ignore_index=True)
                    rows = []
                    
    except Exception as e:
        logging.error(f"Error processing key {key}: {str(e)}")
        continue

# Append any remaining rows
if rows:
    temp_df = pd.DataFrame(rows)
    dataframe_go_panther_complete_slim = pd.concat([dataframe_go_panther_complete_slim, temp_df], ignore_index=True)



    
#           dataframe_go_panther_complete_slim -> 
#           separar em linhas com cada id Uniprot -> 
#           join com cada id UNIPROT do arquivo -> {
#                                            path_origin_kog,
#                                            path_origin_cog,
#                                            path_origin_costa
#                                                 }
#           df's da interseção filtrados por proteina ->
#           separar 
#

#separar para cada uniprot
def collapse_values(series):
    vals = list(map(str, series))
    uniq = list(dict.fromkeys(vals))  # remove duplicados mantendo a ordem
    if len(uniq) == 1:
        return uniq[0]
    else:
        return ",".join(uniq)
    
df_panther_exploded = dataframe_go_panther_complete_slim.assign(Uniprot=dataframe_go_panther_complete_slim['Uniprot'].str.split(',')) \
                .explode('Uniprot') \
                .reset_index(drop=True)

#%%

df_putative_lista_ids_df_filtrado_costa_dme = df_origin_sem_putative_costa[df_origin_sem_putative_costa["Code_Gene_DEG"].isin(lista_ids_df_filtrado_costa_dme)]
df_putative_lista_ids_df_filtrado_costa_sce = df_origin_sem_putative_costa[df_origin_sem_putative_costa["Code_Gene_DEG"].isin(lista_ids_df_filtrado_costa_sce)]
df_putative_lista_ids_df_filtrado_kog_cel = df_origin_sem_putative_kog[df_origin_sem_putative_kog["NCBI_REFSEQ_ACC_VERSION"].isin(lista_ids_df_filtrado_kog_cel)]
df_putative_lista_ids_df_filtrado_kog_dme = df_origin_sem_putative_kog[df_origin_sem_putative_kog["NCBI_REFSEQ_ACC_VERSION"].isin(lista_ids_df_filtrado_kog_dme)]
df_putative_lista_ids_df_filtrado_kog_sce = df_origin_sem_putative_kog[df_origin_sem_putative_kog["NCBI_REFSEQ_ACC_VERSION"].isin(lista_ids_df_filtrado_kog_sce)]
df_putative_lista_ids_df_filtrado_costa_total = df_origin_sem_putative_costa[df_origin_sem_putative_costa["Code_Gene_DEG"].isin(lista_ids_df_filtrado_costa_total)]
df_putative_lista_ids_df_filtrado_kog_total = df_origin_sem_putative_kog[df_origin_sem_putative_kog["NCBI_REFSEQ_ACC_VERSION"].isin(lista_ids_df_filtrado_kog_total)]



df_putative_lista_ids_df_filtrado_kog_cel = df_putative_lista_ids_df_filtrado_kog_cel.rename(columns={"UNIPROTKB_PROTEIN_ACC": "Uniprot"})
df_putative_lista_ids_df_filtrado_kog_dme = df_putative_lista_ids_df_filtrado_kog_dme.rename(columns={"UNIPROTKB_PROTEIN_ACC": "Uniprot"})
df_putative_lista_ids_df_filtrado_kog_sce = df_putative_lista_ids_df_filtrado_kog_sce.rename(columns={"UNIPROTKB_PROTEIN_ACC": "Uniprot"})
df_putative_lista_ids_df_filtrado_kog_total = df_putative_lista_ids_df_filtrado_kog_total.rename(columns={"UNIPROTKB_PROTEIN_ACC": "Uniprot"})


df_putative_lista_ids_df_filtrado_costa_dme_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_costa_dme,df_panther_exploded,on="Uniprot",how="left")
df_putative_lista_ids_df_filtrado_costa_sce_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_costa_sce,df_panther_exploded,on="Uniprot",how="left")
df_putative_lista_ids_df_filtrado_kog_cel_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_kog_cel,df_panther_exploded,on="Uniprot",how="left")
df_putative_lista_ids_df_filtrado_kog_dme_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_kog_dme,df_panther_exploded,on="Uniprot",how="left")
df_putative_lista_ids_df_filtrado_kog_sce_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_kog_sce,df_panther_exploded,on="Uniprot",how="left")

#Rever para os pares
df_putative_lista_ids_df_filtrado_costa_total_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_costa_total,df_panther_exploded,on="Uniprot",how="left") 
df_putative_lista_ids_df_filtrado_kog_total_merged_panther = pd.merge(df_putative_lista_ids_df_filtrado_kog_total,df_panther_exploded,on="Uniprot",how="left") 


#%%

# --- Transformação ---

#df_filtrado_costa_dme
#df_filtrado_costa_sce
#df_filtrado_costa_total
#df_filtrado_kog_cel
#df_filtrado_kog_dme
#df_filtrado_kog_sce
#df_filtrado_kog_total

df_filtrado_costa_dme["costa_droso_pos_go"] = df_filtrado_costa_dme["costa_droso_pos_go"].str.split(", ")
df_explodido_costa_dme = df_filtrado_costa_dme.explode("costa_droso_pos_go").reset_index(drop=True)

df_filtrado_costa_sce["costa_cerevisiae_pos_go"] = df_filtrado_costa_sce["costa_cerevisiae_pos_go"].str.split(", ")
df_explodido_costa_sce = df_filtrado_costa_sce.explode("costa_cerevisiae_pos_go").reset_index(drop=True)

df_filtrado_kog_cel["kog_elegans_pos_go"] = df_filtrado_kog_cel["kog_elegans_pos_go"].str.split(", ")
df_explodido_kog_cel = df_filtrado_kog_cel.explode("kog_elegans_pos_go").reset_index(drop=True)

df_filtrado_kog_dme["kog_droso_pos_go"] = df_filtrado_kog_dme["kog_droso_pos_go"].str.split(", ")
df_explodido_kog_dme = df_filtrado_kog_dme.explode("kog_droso_pos_go").reset_index(drop=True)

df_filtrado_kog_sce["kog_cerevisiae_pos_go"] = df_filtrado_kog_sce["kog_cerevisiae_pos_go"].str.split(", ")
df_explodido_kog_sce = df_filtrado_kog_sce.explode("kog_cerevisiae_pos_go").reset_index(drop=True)

#%%

#####calcular total apos os singles
df_filtrado_kog_total["kog_cerevisiae_pos_go"] = df_filtrado_kog_total["kog_cerevisiae_pos_go"].str.split(", ")
df_explodido_kog_total = df_filtrado_kog_total.explode("kog_cerevisiae_pos_go").reset_index(drop=True)

df_filtrado_costa_total["costa_cerevisiae_pos_go"] = df_filtrado_costa_total["costa_cerevisiae_pos_go"].str.split(", ")
df_explodido_costa_total = df_filtrado_costa_total.explode("costa_cerevisiae_pos_go").reset_index(drop=True)


#merge resultado acima + dataframe com a lista sem putatives e já adicionado GO panther
#df_putative_lista_ids_df_filtrado_costa_dme_merged_panther 
#df_putative_lista_ids_df_filtrado_costa_sce_merged_panther 
#df_putative_lista_ids_df_filtrado_kog_cel_merged_panther
#df_putative_lista_ids_df_filtrado_kog_dme_merged_panther
#df_putative_lista_ids_df_filtrado_kog_sce_merged_panther

#%%
#Renomeando coluna
df_explodido_costa_dme = df_explodido_costa_dme.rename(columns={"costa_droso_pos_go": "Code_Gene_DEG"})
df_explodido_costa_sce = df_explodido_costa_sce.rename(columns={"costa_cerevisiae_pos_go": "Code_Gene_DEG"})
df_explodido_kog_cel = df_explodido_kog_cel.rename(columns={"kog_elegans_pos_go": "NCBI_REFSEQ_ACC_VERSION"})
df_explodido_kog_dme = df_explodido_kog_dme.rename(columns={"kog_droso_pos_go": "NCBI_REFSEQ_ACC_VERSION"})
df_explodido_kog_sce = df_explodido_kog_sce.rename(columns={"kog_cerevisiae_pos_go": "NCBI_REFSEQ_ACC_VERSION"})

df_explodido_costa_total = df_explodido_costa_total.rename(columns={"costa_cerevisiae_pos_go": "Code_Gene_DEG"})
df_explodido_kog_total = df_explodido_kog_total.rename(columns={"kog_cerevisiae_pos_go": "NCBI_REFSEQ_ACC_VERSION"})


dataframe_final_1 = pd.merge(df_putative_lista_ids_df_filtrado_costa_dme_merged_panther, df_explodido_costa_dme,on="Code_Gene_DEG",how="inner")
dataframe_final_2 = pd.merge(df_putative_lista_ids_df_filtrado_costa_sce_merged_panther, df_explodido_costa_sce,on="Code_Gene_DEG",how="inner")
dataframe_final_3 = pd.merge(df_putative_lista_ids_df_filtrado_kog_cel_merged_panther, df_explodido_kog_cel,on="NCBI_REFSEQ_ACC_VERSION",how="inner")
dataframe_final_4 = pd.merge(df_putative_lista_ids_df_filtrado_kog_dme_merged_panther, df_explodido_kog_dme,on="NCBI_REFSEQ_ACC_VERSION",how="inner")
dataframe_final_5 = pd.merge(df_putative_lista_ids_df_filtrado_kog_sce_merged_panther, df_explodido_kog_sce,on="NCBI_REFSEQ_ACC_VERSION",how="inner")

#fazer para todos
dataframe_final_6 = pd.merge(df_putative_lista_ids_df_filtrado_costa_total_merged_panther, df_explodido_costa_total,on="Code_Gene_DEG",how="inner")
dataframe_final_7 = pd.merge(df_putative_lista_ids_df_filtrado_kog_total_merged_panther, df_explodido_kog_total,on="NCBI_REFSEQ_ACC_VERSION",how="inner")


#CONTAGEM FINAL DAS PROTEINAS DO MANSONI

proteinas_mansoni_costa_dme =  set(id_ for sublist in dataframe_final_1["protein_mansoni"].str.split(", ") for id_ in sublist)
proteinas_mansoni_costa_sce =  set(id_ for sublist in dataframe_final_2["protein_mansoni"].str.split(", ") for id_ in sublist)
proteinas_mansoni_kog_cel =  set(id_ for sublist in dataframe_final_3["protein_mansoni"].str.split(", ") for id_ in sublist)
proteinas_mansoni_kog_dme =  set(id_ for sublist in dataframe_final_4["protein_mansoni"].str.split(", ") for id_ in sublist)
proteinas_mansoni_kog_sce =  set(id_ for sublist in dataframe_final_5["protein_mansoni"].str.split(", ") for id_ in sublist)
proteinas_mansoni_costa_total =  set(id_ for sublist in dataframe_final_6["protein_mansoni"].str.split(", ") for id_ in sublist)
proteinas_mansoni_kog_total =  set(id_ for sublist in dataframe_final_7["protein_mansoni"].str.split(", ") for id_ in sublist)


#Colapsando por Orthogroup(Orthofinder)
#Analise sobre GO
def collapse_unique(x):
    # Se todos valores forem iguais, retorna o único valor
    if x.nunique() == 1:
        return x.iloc[0]
    # Caso contrário, junta todos (em string, mas pode ser lista se preferir)
    return ",".join(map(str, x))


#df_colapsado_costa_dme = dataframe_final_1.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))
#df_colapsado_costa_sce = dataframe_final_2.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))
#df_colapsado_costa_total = dataframe_final_6.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))
#df_colapsado_kog_cel = dataframe_final_3.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))
#df_colapsado_kog_dme = dataframe_final_4.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))
#df_colapsado_kog_sce = dataframe_final_5.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))
#df_colapsado_kog_total = dataframe_final_7.groupby("Orthogroup", sort=False).agg(lambda x: ",".join(map(str, x)))

#%%
df_colapsado_costa_dme = dataframe_final_1.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)
df_colapsado_costa_sce = dataframe_final_2.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)
df_colapsado_costa_total = dataframe_final_6.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)
df_colapsado_kog_cel = dataframe_final_3.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)
df_colapsado_kog_dme = dataframe_final_4.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)
df_colapsado_kog_sce = dataframe_final_5.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)
df_colapsado_kog_total = dataframe_final_7.groupby("Orthogroup", sort=False, as_index=False).agg(collapse_unique)


#%%
df_colapsado_kog_sce['id_count'] = df_colapsado_kog_sce['NCBI_REFSEQ_ACC_VERSION'].str.count(',') + 1
df_colapsado_kog_sce['seq_count'] = df_colapsado_kog_sce['NCBI_SEQUENCE'].str.count(',') + 1

# Crie a máscara booleana para identificar as linhas com erro
mask_erros = df_colapsado_kog_sce['id_count'] != df_colapsado_kog_sce['seq_count']

# Use a máscara booleana diretamente para filtrar o DataFrame
# O operador '~' inverte a máscara, mantendo apenas as linhas corretas
df_colapsado_kog = df_colapsado_kog_sce[~mask_erros].copy()

# Continue com o restante do seu código
df_final = pd.DataFrame()
df_final['organism_id'] = df_colapsado_kog['NCBI_REFSEQ_ACC_VERSION'].apply(lambda x: x.split(','))
df_final['sequence'] = df_colapsado_kog['NCBI_SEQUENCE'].apply(lambda x: x.split(','))

df_expanded = df_final.explode(['organism_id', 'sequence']).reset_index(drop=True)

df_expanded['organism_id'] = df_expanded['organism_id'].str.strip()
df_expanded['sequence'] = df_expanded['sequence'].str.strip()

fasta_file = 'kog_sce_final.fasta'

with open(fasta_file, 'w') as f:
    for index, row in df_expanded.iterrows():
        header = f">{row['organism_id']}\n"
        sequence = f"{row['sequence']}\n"
        f.write(header)
        f.write(sequence)

print(f"Arquivo '{fasta_file}' criado com sucesso!")
#%%

#df's limpos derivados 

df_derivado_costa_dme = df_colapsado_costa_dme[["Orthogroup","Code_Gene_DEG","Gene","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


df_derivado_costa_sce = df_colapsado_costa_sce[["Orthogroup","Code_Gene_DEG","Gene","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


df_derivado_costa_total = df_colapsado_costa_total[["Orthogroup","Code_Gene_DEG","Gene","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


df_derivado_kog_cel = df_colapsado_kog_cel[["Orthogroup","KOG","KOG_OL.FUNC","KOG_DESCRICAO","NCBI_REFSEQ_ACC_VERSION","NCBI_GENE","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


df_derivado_kog_dme =  df_colapsado_kog_dme[["Orthogroup","KOG","KOG_OL.FUNC","KOG_DESCRICAO","NCBI_REFSEQ_ACC_VERSION","NCBI_GENE","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


df_derivado_kog_sce =  df_colapsado_kog_sce[["Orthogroup","KOG","KOG_OL.FUNC","KOG_DESCRICAO","NCBI_REFSEQ_ACC_VERSION","NCBI_GENE","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


df_derivado_kog_total =  df_colapsado_kog_total[["Orthogroup","KOG","KOG_OL.FUNC","KOG_DESCRICAO","NCBI_REFSEQ_ACC_VERSION","NCBI_GENE","Uniprot","GO_ID","GO_SYMBOL","GO_QUALIFIER",
                                                "GO_TERM","GO_NAME","GO_PROTEIN_NAME","Panther_GO_id","Term_Label","protein_mansoni"]]


#%%
proteinas_mansoni_cog = set(df_filtrado_OU["protein_mansoni"])

#Lista Final
set_final = proteinas_mansoni_costa_total&\
            proteinas_mansoni_kog_total
            


            #proteinas_mansoni_costa_dme&\
            #proteinas_mansoni_costa_sce&\
            #proteinas_mansoni_kog_sce&\
            #proteinas_mansoni_kog_dme&\
            #proteinas_mansoni_kog_sce

print(len(set_final))
print(set_final)
#%%
# Analise exploratoria de dados

#costa
def contar_ids(serie):
    return (
        serie.dropna()                              
        .str.split(",")                             
        .explode()        
        .loc[lambda s: s.str.lower() != "nan"]      
        .value_counts()                             
    )

#%%
dataframes = {
    "df_costa_dme": df_derivado_costa_dme,
    "df_costa_sce": df_derivado_costa_sce,
    "df_costa_total": df_derivado_costa_total,
    "df_kog_cel": df_derivado_kog_cel,
    "df_kog_dme": df_derivado_kog_dme,
    "df_kog_sce": df_derivado_kog_sce,
    "df_kog_total": df_derivado_kog_total,
}

contagem_go_id = {}
#%%
for nome, df in dataframes.items():


    go_counts = contar_ids(df["GO_ID"])
    panther_counts = contar_ids(df["Panther_GO_id"])

    contagem_go_id[nome] = {
        "GO_ID": go_counts,
        "Panther_GO_id": panther_counts,
        "GO_ID_MAX": go_counts.sum(),
        "Panther_GO_id_MAX": panther_counts.sum(),
        "GO_ID_MEAN":go_counts.mean(),
        "Panther_GO_id_MEAN":panther_counts.mean(),
        "GO_ID_std":go_counts.std(),
        "Panther_GO_id_std":panther_counts.std(),
    }

#Contabilizar as incidências de proteínas que aparecem em todos 
# os dataframes de organismos processados após a interseção 




#passar outro orthofinder?



'''

special_rules = {
    "Origem_organismo": collapse_values,
    "Organismo": collapse_values,       
    "Statistical_Test_Type": collapse_values,
    "Data_Version": collapse_values,
    "Plus_Minus": collapse_values,          
}


# Função default
default_rule = lambda x: ",".join(map(str, x))

# Construir agg_rules completo
agg_rules = {col: special_rules.get(col, default_rule) for col in df_panther_exploded.columns if col != "Uniprot"}

# Aplicar groupby com regras
df_grouped = df_panther_exploded.groupby("Uniprot", as_index=False).agg(agg_rules)


df_origin_sem_putative_cog = df_origin_sem_putative_cog.rename(columns={"UNIPROTKB_PROTEIN_ACC":"Uniprot"})
df_origin_sem_putative_kog = df_origin_sem_putative_kog.rename(columns={"UNIPROTKB_PROTEIN_ACC":"Uniprot"})

df_costa_merged = df_origin_sem_putative_costa.merge(df_grouped, on="Uniprot",how="left")
df_cog_merged = df_origin_sem_putative_cog.merge(df_grouped, on="Uniprot",how="left")
df_kog_merged = df_origin_sem_putative_kog.merge(df_grouped, on="Uniprot",how="left")

#Filtrados da interseção

df_filter_costa_dme_exploded = df_filtrado_costa_dme.assign(costa_droso_pos_go=df_filtrado_costa_dme['costa_droso_pos_go'].str.split(', ')) \
                .explode('costa_droso_pos_go') \
                .reset_index(drop=True)
                
df_filter_costa_sce_exploded = df_filtrado_costa_sce.assign(costa_cerevisiae_pos_go = df_filtrado_costa_sce['costa_cerevisiae_pos_go'].str.split(', ')) \
                .explode('costa_cerevisiae_pos_go') \
                .reset_index(drop=True)  
                
#df_filter_costa_total_exploded = df_filtrado_costa_total.assign(Uniprot=df_filtrado_costa_total['costa_cerevisiae_pos_go'].str.split(',')) \
#                .explode('costa_cerevisiae_pos_go') \
#                .reset_index(drop=True)                 

df_filter_kog_cel_exploded = df_filtrado_kog_cel.assign(kog_elegans_pos_go=df_filtrado_kog_cel['kog_elegans_pos_go'].str.split(', ')) \
                .explode('kog_elegans_pos_go') \
                .reset_index(drop=True)  
df_filter_kog_dme_exploded = df_filtrado_kog_dme.assign(kog_droso_pos_go=df_filtrado_kog_dme['kog_droso_pos_go'].str.split(', ')) \
                .explode('kog_droso_pos_go') \
                .reset_index(drop=True)  
df_filter_kog_sce_exploded = df_filtrado_kog_sce.assign(kog_cerevisiae_pos_go=df_filtrado_kog_sce['kog_cerevisiae_pos_go'].str.split(', ')) \
                .explode('kog_cerevisiae_pos_go') \
                .reset_index(drop=True)   

#renome para o join
#rever esta parte 

df_filter_costa_dme_exploded = df_filter_costa_dme_exploded.rename(columns={"costa_droso_pos_go":"Code_Organism"})
df_filter_costa_sce_exploded = df_filter_costa_sce_exploded.rename(columns={"costa_cerevisiae_pos_go":"Code_Organism"})
df_filter_kog_cel_exploded = df_filter_kog_cel_exploded.rename(columns={"kog_elegans_pos_go":"NCBI_REFSEQ_ACC_VERSION"})
df_filter_kog_dme_exploded = df_filter_kog_dme_exploded.rename(columns={"kog_droso_pos_go":"NCBI_REFSEQ_ACC_VERSION"})
df_filter_kog_sce_exploded = df_filter_kog_sce_exploded.rename(columns={"kog_cerevisiae_pos_go":"NCBI_REFSEQ_ACC_VERSION"})


df_final_merged_intersecao_sem_putative_costa_1   = df_costa_merged.merge(df_filter_costa_dme_exploded, on="Code_Organism",how="inner")
df_final_merged_intersecao_sem_putative_costa_2   = df_costa_merged.merge(df_filter_costa_sce_exploded, on="Code_Organism",how="inner")
                                                
df_final_merged_intersecao_sem_putative_kog_1 = df_kog_merged.merge(df_filter_kog_cel_exploded,on="NCBI_REFSEQ_ACC_VERSION",how="inner")
df_final_merged_intersecao_sem_putative_kog_2 = df_kog_merged.merge(df_filter_kog_dme_exploded,on="NCBI_REFSEQ_ACC_VERSION",how="inner")
df_final_merged_intersecao_sem_putative_kog_3 = df_kog_merged.merge(df_filter_kog_sce_exploded,on="NCBI_REFSEQ_ACC_VERSION",how="inner")



#separar as proteínas do mansoni

# kog: cel -> 732; dme -> 546; sce -> 550;
# costa: dme -> 0; sce -> 0;

lista_ids_kog_1 = df_final_merged_intersecao_sem_putative_kog_1["protein_mansoni"].str.split(", ").explode().tolist()
lista_ids_kog_2 = df_final_merged_intersecao_sem_putative_kog_2["protein_mansoni"].str.split(", ").explode().tolist()
lista_ids_kog_3 = df_final_merged_intersecao_sem_putative_kog_3["protein_mansoni"].str.split(", ").explode().tolist()






'''







'''
    else:
        origem = name_temp[0].split('_')[0]
        organismo = value["overrepresentation"]["upload_lists"]["input_list"]["organism"]
        Annotation_Type =  value["overrepresentation"]["annotation_type"]
        Bonferroni_Correction_value = value["overrepresentation"]["bonferroni_correction"]
        Statistical_Test_Type = value["overrepresentation"]["test_type"]
        Data_Version = value["overrepresentation"]["data_version_release_date"]
        
        for grupos in value["overrepresentation"]["group"]:
            for grupo in grupos["result"]:
                for result_dict in grupo:
                    dataframe_go_panther_complete_slim.append(
                        {
                            "Origem_organismo":origem,
                            "Organismo":organismo,
                            "Annotation_Type":Annotation_Type,
                            "Statistical_Test_Type":Statistical_Test_Type,
                            "Bonferroni_Correction_value":Bonferroni_Correction_value,
                            "Data_Version":Data_Version,
                            "Panther_GO_id":result_dict["term"]["id"],
                            "Term_Label":result_dict["term"]["label"],
                            "Term_Level":result_dict["term"]["level"],
                            "QTD_proteins_in_list":result_dict["input_list"]["number_in_list"],
                            "Expected":result_dict["input_list"]["expected"],
                            "Fold_Enrichment":result_dict["input_list"]["fold_enrichment"],
                            "pValue":result_dict["input_list"]["pValue"],
                            "Plus_Minus":result_dict["input_list"]["plus_minus"],
                            "QTD_proteins_in_Reference":result_dict["number_in_reference"],
                            "Protein_IDs":",".join(result_dict["input_list"]["mapped_id_list"]["mapped_id"])#result_dict["input_list"]["mapped_id_list"]["mapped_id"],
                            }, ignore_index=True)
'''                    

'''  
# Initialize the directed graph
for arvore in dicionario_analises_arquivos_panther.keys():
  if arvore.endswith("go_complete.json"):
      G = nx.DiGraph()
      root = "overrepresentation"
      G.add_node(root, level=-1)

      # Get the "group" data
      group_data = dicionario_analises_arquivos_panther[arvore].get("overrepresentation", {}).get("group", {})

      # Ensure group_data is a list; if it's a dictionary, wrap it in a list
      if isinstance(group_data, dict):
          group_data = [group_data]
      elif not isinstance(group_data, list):
          print(f"Error: 'group' is not a list or dict in {arvore}, got {type(group_data)}: {group_data}")
          group_data = []

      # Iterate through groups
      for group in group_data:
          
          # Get the "result" data
          result_data = group.get("result", [])
          
          # Ensure result_data is a list
          if isinstance(result_data, str):
              print(f"Error: 'result' is a string in {arvore}, got: {result_data}")
              continue
          elif not isinstance(result_data, list):
              print(f"Error: 'result' is not a list in {arvore}, got {type(result_data)}: {result_data}")
              continue
          
          # Iterate through results
          for result in result_data:
              # Ensure "term" exists and is a dictionary
              term = result.get("term", {})
              if not isinstance(term, dict):
                  print(f"Error: 'term' is not a dictionary in {arvore}, got {type(term)}: {term}")
                  continue
                  
              # Extract term details with defaults
              term_id = term.get("id", "Unknown_ID")
              term_label = term.get("label", "Unknown_Label")
              level = term.get("level", -1)
              label = f"{term_id} - {term_label}"

              # Add node
              G.add_node(label, level=level)

              # Connect by hierarchy
              if level == 0:
                  G.add_edge(root, label)
              else:
                  # Find the parent node (immediate previous level)
                  parent = None
                  for n, attrs in G.nodes(data=True):
                      if attrs.get("level") == level - 1:
                          parent = n
                          break
                  if parent:
                      G.add_edge(parent, label)

      # Layout and plot
      try:
          pos = graphviz_layout(G, prog="dot")
      except ImportError:
          print("Warning: graphviz_layout not available, falling back to spring_layout")
          pos = nx.spring_layout(G, seed=33,k=0.5)

      # Plot the graph
      plt.figure(figsize=(16, 12))
      nx.draw(
          G, pos,
          with_labels=True,
          node_size=4500,
          node_color="lightblue",
          font_size=8,
          font_weight="bold",
          arrows=True
      )
      # Título
      plt.title(f"GO Term Hierarchy para {arvore}", pad=20)

      # Ajustar margens automaticamente
      plt.tight_layout()
      plt.margins(x=0.1, y=0.1)

      # Mostrar o gráfico
      plt.show()
'''
