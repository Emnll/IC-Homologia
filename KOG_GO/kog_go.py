# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 11:42:30 2025

@author: Jean
"""

import pandas as pd
import os 
import requests, sys
import json
import time
import datetime

def fetchQuickGo(uniprot_id):
#Função que realiza o envio de requisição para busca de informações de Ontologia de Gene utilizando a API QuickGo.
#Aceita como parâmetro uma string de ids separados por vírgula e retorna as informações acerca da Ontologia, sendo elas:
#    GO_ID => Identificador do Gene Ontology,
#    GO_SYMBOL => Correspondência do Gene Product ID ,
#    "GO_QUALIFIER" => Define a relação entre Gene Product e GO term,
#    "GO_TERM"=> Identificador único da função ,
#    "GO_EVIDENCECODE" => É a representação da atibuição dos GO terms as proteínas com base em diferentes evidências.
#Há a requisição da string porém o body da requisição é dividido em páginas pela API, requerindo uma iteração das páginas (novas requisições modificando as páginas)
    pagina = "1"
    requestURL_page1 = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?includeFields=goName,taxonName,name&geneProductId=" + uniprot_id + "&limit=200" +"&page="+ pagina
    
    r_pg1 = requests.get(requestURL_page1, headers={ "Accept" : "application/json"})

    if not r_pg1.ok:
      r_pg1.raise_for_status()
      sys.exit()
    
    #Dado o body da resposta, retira os resultados da primeira página com intuito de formatação e adiciona posteriormente 
    result_pg1 = r_pg1.json()["results"]
    responseBody_pg = (r_pg1.json())
    responseBody_pg.pop("results")   
    responseBody_pg["results"] = []
    responseBody_pg["results"].append(result_pg1)

    if responseBody_pg['pageInfo']['total'] > 1:
        
        for i in range(2,responseBody_pg['pageInfo']['total'] + 1):
            time.sleep(1)
            pagina = str(i)
            requestURL_pagen = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?includeFields=goName,taxonName,name&geneProductId=" + uniprot_id + "&limit=200" +"&page="+ pagina
            
            r_pgn = requests.get(requestURL_pagen, headers={ "Accept" : "application/json"})

            if not r_pgn.ok:
              r_pgn.raise_for_status()
              sys.exit()
              
            pg_result = r_pgn.json()
            responseBody_pg["results"].append(pg_result["results"])
    
    return responseBody_pg

inicio = datetime.datetime.now()
#Leitura do caminho 
diretorio_atual = os.getcwd()
resultado_arq = os.path.join(diretorio_atual,"resultado_csv.csv")


df_kog =  pd.read_csv(resultado_arq,
            sep=",")

#Criação do dataframe com retirada de proteínas não indexadas com o arquivo collab do NCBI x UniProtKB.
#A retirada destas proteínas são necessárias pois mesmo que elas já estejam atualizadas os registros das mesmas
#foram suprimidas pelo processamento padrão de anotação de genomas do ncbi. 

df_kog_drop_uniprotkb = df_kog[df_kog["UNIPROTKB_PROTEIN_ACC"].notna()]

uniprots = df_kog_drop_uniprotkb.loc[:,"UNIPROTKB_PROTEIN_ACC"].tolist()
#Preparo das strings para utilização da API QuickGo
#A execução fica entorno de 01:10hr 
batches_uniprots = []
n = 200
string = ""
body_pages = {}

for i in range(0,len(uniprots),n):
    batch_string = uniprots[i:i+n]
    string = ','.join(map(str,batch_string))
    batches_uniprots.append(string)
    
for batch in range(0,len(batches_uniprots)):
    batch_results = fetchQuickGo(batches_uniprots[batch])
    body_pages[batch] = batch_results['results']

final_quickgo = datetime.datetime.now()
diff_quickgo = final_quickgo - inicio
total_segundos_quickgo = diff_quickgo.total_seconds()
horas = int(total_segundos_quickgo// 3600)
minutos = int((total_segundos_quickgo % 3600) // 60)
segundos = int(total_segundos_quickgo% 60)
fractions = int((total_segundos_quickgo - int(total_segundos_quickgo)) * 1000)
data_hora_formatada = f"{horas:02}:{minutos:02}:{segundos:02}:{fractions:03}"
print("#####################################################\n")
print(f"Tempo total das requisições QuickGO: {data_hora_formatada}")
#criação novo dataframe para futuro merge
# POSSIVEL FUTURO TRABALHO É VER POR INFORMAÇÃO SEMANTICA PELO GO REFERENCE
dataframe_GO = pd.DataFrame(columns=[
    "UNIPROTKB_PROTEIN_ACC",
    "GO_ID",
    "GO_SYMBOL",
    "GO_QUALIFIER",
    "GO_TERM",
    "GO_EVIDENCECODE"
    ])

lista_go_final = []

for batch in body_pages.values():
    for lista_page in batch:
        for dicionario_page in lista_page:
            lista_go_final.append({
                "UNIPROTKB_PROTEIN_ACC":dicionario_page["geneProductId"].split(':')[1],
                "GO_ID":dicionario_page["goId"],
                "GO_SYMBOL":dicionario_page["symbol"],
                "GO_QUALIFIER":dicionario_page["qualifier"],
                "GO_TERM":dicionario_page["goAspect"],
                "GO_EVIDENCECODE":dicionario_page["evidenceCode"],
                "GO_NAME": dicionario_page["goName"],
                "GO_PROTEIN_NAME": dicionario_page["name"]
                })


#criar um dicionario  para organizar mesmas entradas de uniprots
dict_org_go = {}

for i in range(0,len(lista_go_final)):
    dict_inList = lista_go_final[i]
    uniprotkb = dict_inList["UNIPROTKB_PROTEIN_ACC"]
    if uniprotkb in dict_org_go:
        # ja existe: concatena cada campo com virgula
        dict_org_go[uniprotkb]["GO_ID"] += ',' + dict_inList["GO_ID"]
        dict_org_go[uniprotkb]["GO_QUALIFIER"] += ',' + dict_inList["GO_QUALIFIER"]
        dict_org_go[uniprotkb]["GO_TERM"] += ',' + dict_inList["GO_TERM"]
        dict_org_go[uniprotkb]["GO_EVIDENCECODE"] += ',' + dict_inList["GO_EVIDENCECODE"]
        dict_org_go[uniprotkb]["GO_NAME"] += ',' + dict_inList["GO_NAME"]
    
    else:
        # primeira vez: cria o dicionario com os valores iniciais
        dict_org_go[uniprotkb] = {
            "GO_ID": dict_inList["GO_ID"],
            "GO_SYMBOL": dict_inList["GO_SYMBOL"],
            "GO_QUALIFIER": dict_inList["GO_QUALIFIER"],
            "GO_TERM": dict_inList["GO_TERM"],
            "GO_EVIDENCECODE": dict_inList["GO_EVIDENCECODE"],
            "GO_NAME": dict_inList["GO_NAME"],
            "GO_PROTEIN_NAME": dict_inList["GO_PROTEIN_NAME"] 
        }


list_formatacao_go = []

for k, v in dict_org_go.items():
    list_formatacao_go.append({
        "UNIPROTKB_PROTEIN_ACC":k,
        "GO_ID":v["GO_ID"],
        "GO_SYMBOL":v["GO_SYMBOL"],
        "GO_QUALIFIER":v["GO_QUALIFIER"],
        "GO_TERM":v["GO_TERM"],
        "GO_EVIDENCECODE":v["GO_EVIDENCECODE"],
        "GO_NAME": v["GO_NAME"],
        "GO_PROTEIN_NAME": v["GO_PROTEIN_NAME"]
        })


#Dado um UniprotKB, a relação entre GO_ID,GO_QUALIFIER,GO_TERM,GO_EVIDENCECODE,GO_NAME segue a ordem encontrada nas strigs separadas por ",".
#Ex.: 
#    O18165:
#    
#      GO_ID: 
#        GO:0004198         / GO:0006508
#      GO_SYMBOL:
#        clpr-1
#      GO_QUALIFIER:
#        enables            / involved_in
#      GO_TERM:
#        molecular_function / biological_process
#      GO_EVIDENCECODE:
#        ECO:0000256        / ECO:0000256
#      GO_NAME:
#        calcium-dependent cysteine-type endopeptidase activity        / proteolysis
#      GO_PROTEIN_NAME:
#        Calpain catalytic domain-containing protein
        
    
dataframe_temp = pd.DataFrame(list_formatacao_go)

#Merge dataframe que retirou-se os dados suprimidos e o dataframe formatado do GO

dataframe_left_GO = pd.merge(df_kog_drop_uniprotkb,dataframe_temp, on="UNIPROTKB_PROTEIN_ACC", how='left')
dataframe_right_GO = pd.merge(df_kog_drop_uniprotkb,dataframe_temp, on="UNIPROTKB_PROTEIN_ACC", how='right')


#Escrita dos arquivos resultates:
#   "dataframe_temp", que detém informações separadas de GO sobre UNIPROTKD id do grupo de KOG e; 
#
#   "dataframe_left_GO": 22470, dataframe final com todas as informações do merge entre tabela GO(right) e tabela kog(left) com proteínas suprimidas retiradas;
#
#   "dataframe_right_GO": 15117, dataframe final com todas as informações onde o merge da última tabela
#         é dado pela tabela com informação do GO(left) e tabela kog(right). 

nome_arquivo1 = os.path.join(diretorio_atual,"kog_go_informations.csv")
nome_arquivo2 = os.path.join(diretorio_atual,"Resultado_kog_GO_left.csv")
nome_arquivo3 = os.path.join(diretorio_atual,"Resultado_kog_GO_right.csv")

dataframe_temp.to_csv(nome_arquivo1,index=False)
dataframe_left_GO.to_csv(nome_arquivo2,index=False)
dataframe_right_GO.to_csv(nome_arquivo3,index=False)

final = datetime.datetime.now()
diff = final - inicio
total_segundos = diff.total_seconds()
horas = int(total_segundos // 3600)
minutos = int((total_segundos % 3600) // 60)
segundos = int(total_segundos % 60)
fractions = int((total_segundos - int(total_segundos)) * 1000)
data_hora_formatada = f"{horas:02}:{minutos:02}:{segundos:02}:{fractions:03}"
print("#####################################################\n")
print(f"Tempo total de execução: {data_hora_formatada}")