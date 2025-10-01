# -*- coding: utf-8 -*-
"""
Created on Wed May  7 22:54:29 2025

@author: Jean
"""

import os
import pandas as pd
import requests, sys, json
import re
import time
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry
import datetime

inicio = datetime.datetime.now()
#base_url = "https://rest.uniprot.org/"
POLLING_INTERVAL = 3
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


API_URL = "https://rest.uniprot.org"

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


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)


diretorio_atual = os.getcwd()

#txt_uniprot = 'gene_refseq_uniprotkb_collab.txt' 
txt_uniprot = r"C:\Users\Jean\Desktop\periodo 2024.2\IC\apresentação 3\TESTE_CODIGO_IC\gene_refseq_uniprotkb_collab.txt"

#url_cog = r"dataframe_cog_completo_ou.csv"
url_cog = r"C:\Users\Jean\Desktop\periodo 2024.2\IC\Entrega orthofinder\código\RESULTADOS FINAIS\dataframe_cog_completo_ou.csv"

df_csv = pd.read_csv(url_cog, sep=',',low_memory=False)


chunksize = 10**6 #Alterar o valor caso esteja com dificuldades em processar o dado (chunks menores utilizaram menos memória porém processo dura mais)
df_csv['UNIPROTKB_PROTEIN_ACC'] = None

#há retirada de duplicatas onde apenas os valores de footprints que diferenciavam as duas entradas.

df_csv = df_csv.drop_duplicates(subset='Protein_ID', keep='first')

#Troca do padrão de id para comparar com o uniprot

df_csv['Protein_ID'] = df_csv['Protein_ID'].str.replace(r'_(?=[^_]*$)', '.', regex=True)

NCBI_REFSEQ_ACCVER = set(df_csv['Protein_ID'])


path_uniprot = txt_uniprot #os.path.join(diretorio_atual,txt_uniprot)

for chunk in pd.read_csv(
        path_uniprot,
        chunksize=chunksize,
        sep='\t',
        usecols=['#NCBI_protein_accession','UniProtKB_protein_accession'],
        dtype= str
        ):
    chunk_filtro = chunk[chunk['#NCBI_protein_accession'].isin(NCBI_REFSEQ_ACCVER)]
    
    for _,row in chunk_filtro.iterrows():
        ncbi_acc = row['#NCBI_protein_accession']
        uniprot_acc = row['UniProtKB_protein_accession']
        df_csv.loc[df_csv['Protein_ID'] == ncbi_acc,'UNIPROTKB_PROTEIN_ACC'] = uniprot_acc

final_collab = datetime.datetime.now()
diff_collab = final_collab - inicio
total_segundos_collab = diff_collab.total_seconds()
horas = int(total_segundos_collab// 3600)
minutos = int((total_segundos_collab % 3600) // 60)
segundos = int(total_segundos_collab% 60)
fractions = int((total_segundos_collab - int(total_segundos_collab)) * 1000)
data_hora_formatada = f"{horas:02}:{minutos:02}:{segundos:02}:{fractions:03}"
print("#####################################################\n")
print(f"Tempo total criação dataframe Uniprot x NCBI: {data_hora_formatada}")
UNIPROTKB_NONE = df_csv['UNIPROTKB_PROTEIN_ACC'].isnull()

ID_BUSCA = df_csv.loc[UNIPROTKB_NONE,'Protein_ID'].tolist()


job_id = submit_id_mapping(
    from_db="RefSeq_Protein", to_db="UniProtKB", ids= ID_BUSCA
)
if check_id_mapping_results_ready(job_id):
    link = get_id_mapping_results_link(job_id)
    results = get_id_mapping_results_search(link)
    # Equivalently using the stream endpoint which is more demanding
    # on the API and so is less stable:
    # results = get_id_mapping_results_stream(link)

#failedIds => ids extras numa relação 1:n, tendo sido utilizado a melhor correspodência (UniprotKB) para BD (e.g. Uniref90,Uniref50 e etc).
#    lista completa pode ser achada em: https://rest.uniprot.org/configure/uniprotkb/search-fields
#suggestedIds => São atribuídos valores Uniparc para algumas das entradas. Os ids Uniparcs(2101 identificadores não utilizados) não serviram para buscar os GO visto a natureza do identificador
#   serve para rastreio do histórico das mudanças da sequência em diferentes databases
#    O conteudo na íntegra, pode ser acessado em: https://www.uniprot.org/help/uniparc
#obsoleteCout => contador de id obsoletos no Uniprot adicionandos a entrada
#results => objeto visado onde estão as proteínas que não foram achadas pelo arquivo de colaboração mas pela busca via API
#   *Há algumas informações de GO que são retornadas na parte de referência cruzada, porém não será utilizado, sendo necessário ainda a buscca via QuickGO
#        pois há informações faltantes como: GO_symbol,GO_qualifier e etc.

lista_busca_uniprot = []
[results["results"]]
for result in results["results"]:
    lista_busca_uniprot.append({
        "Protein_ID":result["from"],
        "UNIPROTKB_PROTEIN_ACC":result["to"]["primaryAccession"]
        })

pd_busca_uniprot = pd.DataFrame(lista_busca_uniprot)

df_uniprot_final = pd.merge(df_csv,pd_busca_uniprot,on="Protein_ID",how='left')
df_uniprot_final['UNIPROTKB_PROTEIN_ACC'] = df_uniprot_final["UNIPROTKB_PROTEIN_ACC_x"].combine_first(df_uniprot_final['UNIPROTKB_PROTEIN_ACC_y'])
df_uniprot_final.drop(columns=["UNIPROTKB_PROTEIN_ACC_x","UNIPROTKB_PROTEIN_ACC_y"],inplace=True)

#retirada das entradas com Uniprot nulos
df_uniprot_final_dropped = df_uniprot_final[df_uniprot_final["UNIPROTKB_PROTEIN_ACC"].notna()]

uniprots = df_uniprot_final_dropped.loc[:,"UNIPROTKB_PROTEIN_ACC"].tolist()


#Preparo das strings para utilização da API QuickGo
#A execução fica entorno de 01:51 hr
inicio_quickgo = datetime.datetime.now()
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
diff_quickgo = final_quickgo - inicio_quickgo
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

dataframe_left_GO = pd.merge(df_uniprot_final_dropped,dataframe_temp, on="UNIPROTKB_PROTEIN_ACC", how='left')
dataframe_right_GO = pd.merge(df_uniprot_final_dropped,dataframe_temp, on="UNIPROTKB_PROTEIN_ACC", how='right')

nome_arquivo1 = os.path.join(diretorio_atual,"cog_go_informations.csv")
nome_arquivo2 = os.path.join(diretorio_atual,"Resultado_cog_GO_left.csv")
nome_arquivo3 = os.path.join(diretorio_atual,"Resultado_cog_GO_right.csv")

dataframe_temp.to_csv(nome_arquivo1,index=False)
dataframe_left_GO.to_csv(nome_arquivo2,index=False)
dataframe_right_GO.to_csv(nome_arquivo3,index=False)
            
#nome_csv= os.path.join(join_pastafinal, '.csv')
#df_csv.to_csv(nome_csv,header=True,index=False)