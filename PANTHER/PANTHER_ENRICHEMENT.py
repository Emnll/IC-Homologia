import requests
import pandas as pd
import os
import re
import time

# URLs da API do PANTHER
#url_enrich = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep"
#url_validate = "http://pantherdb.org/services/oai/tools/validateIds"

'''
# Função para validar IDs de genes
def validate_gene_ids(gene_list, organism_code):
    gene_list_str = ",".join(gene_list)
    print(gene_list_str)
    params = {
        "ids": gene_list_str,
        "taxonId": organism_code
    }
    try:
        response = requests.get(url_validate, params=params)
        response.raise_for_status()
        result = response.json()
        mapped_ids = [id_info.get("id") for id_info in result.get("mapped", [])]
        unmapped_ids = result.get("unmapped", [])
        print(f"Validação para organismo {organism_code}: {len(mapped_ids)} IDs mapeados, {len(unmapped_ids)} IDs não mapeados")
        if unmapped_ids:
            print(f"IDs não mapeados: {unmapped_ids[:10]}...")
        return mapped_ids
    except Exception as e:
        print(f"Erro ao validar IDs para organismo {organism_code}: {str(e)}")
        return []

# Função para realizar enriquecimento por organismo
def perform_enrichment(df, df_name, organism_code, gene_column, annotdataset):
    gene_list = df[gene_column].dropna().unique().tolist()
    if not gene_list:
        print(f"Nenhum gene encontrado para {organism_code} em {df_name}")
        return None

    # Validar IDs
    valid_genes = validate_gene_ids(gene_list, organism_code)
    if not valid_genes:
        print(f"Nenhum gene válido para {organism_code} em {df_name}")
        return None

    print(f"Genes válidos para {df_name} ({organism_code}): {len(valid_genes)} genes")
    gene_list_str = ",".join(valid_genes)
    
    params = {
        "geneInputList": gene_list_str,
        "organism": organism_code,
        "annotDataSet": annotdataset,
        "enrichmentTestType": "FISHER",
        "correction": "FDR"
    }
    
    try:
        response = requests.get(url_enrich, params=params)
        if response.status_code == 200:
            results = response.json().get("results", {}).get("result", [])
            if not results:
                print(f"Nenhum resultado de enriquecimento para {organism_code} em {df_name} com {annotdataset}")
                return None
            enriched_terms = pd.DataFrame([
                {
                    "GO_ID": res["term"]["id"],
                    "GO_Term": res["term"]["label"],
                    "Genes_in_List": res["number_in_list"],
                    "Genes_in_Reference": res["number_in_reference"],
                    "Fold_Enrichment": res["fold_enrichment"],
                    "P_Value": res["pValue"],
                    "FDR": res["fdr"]
                } for res in results
            ])
            enriched_terms = enriched_terms[enriched_terms['FDR'] < 0.05]
            if enriched_terms.empty:
                print(f"Nenhum termo significativo (FDR < 0.05) para {organism_code} em {df_name} com {annotdataset}")
                return None
            output_file = f"enriched_go_results_{df_name}_{organism_code}_{annotdataset.replace(':', '_')}.csv"
            enriched_terms.to_csv(output_file, index=False)
            print(f"Resultados salvos em {output_file}")
            return enriched_terms
        else:
            print(f"Erro na requisição para {organism_code} em {df_name}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Erro ao processar {organism_code} em {df_name}: {str(e)}")
        return None

'''
# Função que extrai todos os termos entre colchetes []
def extrair_organismos_multiplos(definicao):
    return re.findall(r'\[([^\[\]]+)\]', str(definicao))

# Função para verificar se há algum dos esperados
def tem_organismo_esperado(lista_organismos):
    return any(org in lista_organismos for org in organismos_esperados)

# Caminhos dos arquivos
path_KOG_GO = r"KOG_GO\Resultado_kog_GO_left_sem_putative.csv"
path_COSTA_GO = r"COSTA_GO\Resultado_costa_GO_left_sem_putative.csv"

# Diretório atual
diretorio_atual = os.path.split(os.getcwd())[0]
#diretorio_atual = os.getcwd()

# Caminhos completos
caminho_kog = os.path.join(diretorio_atual, path_KOG_GO)
caminho_costa = os.path.join(diretorio_atual, path_COSTA_GO)

# Carregar os dataframes
pd_kog = pd.read_csv(caminho_kog, sep=',')
pd_costa = pd.read_csv(caminho_costa, sep=',')

# Remover duplicatas e filtrar dados
pd_kog_sem_duplicatas = pd_kog.drop_duplicates(subset="NCBI_REFSEQ_ACC_VERSION", keep='first')
pd_costa_sem_notavailable = pd_costa[pd_costa["Seq_Prot"] != "Notavailablenow."]

# Processar KOG
organismos_esperados = [
    "Caenorhabditis elegans",
    "Homo sapiens",
    "Drosophila melanogaster",
    "Saccharomyces cerevisiae S288C"
]

pd_kog_sem_duplicatas["Organismos_extraidos"] = pd_kog_sem_duplicatas["NCBI_DEFINITION"].apply(extrair_organismos_multiplos)
pd_kog_sem_duplicatas["Tem_organismo_esperado"] = pd_kog_sem_duplicatas["Organismos_extraidos"].apply(tem_organismo_esperado)

df_esperados = pd_kog_sem_duplicatas[pd_kog_sem_duplicatas["Tem_organismo_esperado"] == True]
df_nao_esperados = pd_kog_sem_duplicatas[pd_kog_sem_duplicatas["Tem_organismo_esperado"] == False]

df_kog_human = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Homo sapiens" in x)]
df_kog_droso = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Drosophila melanogaster" in x)]
df_kog_yeast = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Saccharomyces cerevisiae S288C" in x)]
df_kog_elegans = df_esperados[df_esperados["Organismos_extraidos"].apply(lambda x: "Caenorhabditis elegans" in x)]

# Processar COSTA
organismos_unicos = pd_costa_sem_notavailable["Organism"].unique()
dfs_por_organismo = {
    org: pd_costa_sem_notavailable[pd_costa_sem_notavailable["Organism"] == org].copy()
    for org in organismos_unicos
}

df_costa_human = dfs_por_organismo.get("Homo sapiens")
df_costa_droso = dfs_por_organismo.get("Drosophila melanogaster")
df_costa_yeast = dfs_por_organismo.get("Saccharomyces cerevisiae")


nome_arquivo_kog_human = os.path.join(diretorio_atual,"Resultado_kog_human_para_enrichment.csv")
nome_arquivo_kog_droso = os.path.join(diretorio_atual,"Resultado_kog_droso_para_enrichment.csv")
nome_arquivo_kog_yeast = os.path.join(diretorio_atual,"Resultado_kog_yeast_para_enrichment.csv")
nome_arquivo_kog_elegans = os.path.join(diretorio_atual,"Resultado_kog_elegans_para_enrichment.csv")
nome_arquivo_costa_human = os.path.join(diretorio_atual,"Resultado_costa_human_para_enrichment.csv")
nome_arquivo_costa_droso = os.path.join(diretorio_atual,"Resultado_costa_droso_para_enrichment.csv")
nome_arquivo_costa_yeast = os.path.join(diretorio_atual,"Resultado_costa_yeast_para_enrichment.csv")


df_kog_human.to_csv(nome_arquivo_kog_human,index=False)
df_kog_droso.to_csv(nome_arquivo_kog_droso,index=False)
df_kog_yeast.to_csv(nome_arquivo_kog_yeast,index=False)
df_kog_elegans.to_csv(nome_arquivo_kog_elegans,index=False)
df_costa_human.to_csv(nome_arquivo_costa_human,index=False)
df_costa_droso.to_csv(nome_arquivo_costa_droso,index=False)
df_costa_yeast.to_csv(nome_arquivo_costa_yeast,index=False)


 



'''
# Definir anotações
annots = [
    "GO:0003674",  # Molecular Function
    "GO:0008150",  # Biological Process
    "GO:0005575"   # Cellular Component
]

# Dicionários para armazenar resultados
resultado_kog_human = {}
resultado_kog_droso = {}
resultado_kog_yeast = {}
resultado_kog_elegans = {}
resultado_costa_human = {}
resultado_costa_droso = {}
resultado_costa_yeast = {}

# Executar enriquecimento
for annotation in annots:
    print(f"Processando anotação {annotation}...")
    if not df_kog_human.empty:
        resultado_kog_human[annotation] = perform_enrichment(df_kog_human, "KOG_HUMAN", '9606', "UNIPROTKB_PROTEIN_ACC", annotation)
        time.sleep(3)
    if not df_kog_droso.empty:
        resultado_kog_droso[annotation] = perform_enrichment(df_kog_droso, "KOG_DROSO", '7227', "UNIPROTKB_PROTEIN_ACC", annotation)
        time.sleep(3)
    if not df_kog_yeast.empty:
        resultado_kog_yeast[annotation] = perform_enrichment(df_kog_yeast, "KOG_YEAST", '559292', "UNIPROTKB_PROTEIN_ACC", annotation)
        time.sleep(3)
    if not df_kog_elegans.empty:
        resultado_kog_elegans[annotation] = perform_enrichment(df_kog_elegans, "KOG_ELEGANS", '6239', "UNIPROTKB_PROTEIN_ACC", annotation)
        time.sleep(3)
    if df_costa_human is not None and not df_costa_human.empty:
        resultado_costa_human[annotation] = perform_enrichment(df_costa_human, "COSTA_HUMAN", '9606', "Uniprot", annotation)
        time.sleep(3)
    if df_costa_droso is not None and not df_costa_droso.empty:
        resultado_costa_droso[annotation] = perform_enrichment(df_costa_droso, "COSTA_DROSO", '7227', "Uniprot", annotation)
        time.sleep(3)
    if df_costa_yeast is not None and not df_costa_yeast.empty:
        resultado_costa_yeast[annotation] = perform_enrichment(df_costa_yeast, "COSTA_YEAST", '559292', "Uniprot", annotation)
        time.sleep(3)
'''