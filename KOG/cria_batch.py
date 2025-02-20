# -*- coding: utf-8 -*-

import datetime
import pandas as pd
from Bio import Entrez
import os
import json
import shutil
import re
import time

email = input("Entrez email(ou usuário anônimo, simplesmente apertando Enter): ").strip() or None
Entrez.email = email
apikey = input("Entrez API key(Caso não use, apenas três queries por segundo será permitido): ").strip() or None
Entrez.api_key= apikey

#funcoes
def limpeza_kog(protein_name):
    #funcao usada para retirar encontros no kog como em: CE02013_2, Onde a proteina se refere apenas a CE02013,
    #sendo "_" apenas para representar a quantidade de vezes da proteina atribuida a um KOG distinto.
    
    if "_" in protein_name:
        return protein_name.split('_')[0]
    return protein_name
def cria_batch(dict_gi):
    #   Retorna um dicionario, onde a chave é o organismo e o valor, uma lista de 
    #outras listas, contendo strings de noo máximo 199 gi's daquele organismos.
    #
    #   Por uma limitação de uso da API e do método Efetch, foi posto no máximo 
    #199 GenbankID para query
    
    dicionario = {}
    n = 199
    for key,gi_list in dict_gi.items():
        sub_listas = []
        for i in range(0,len(gi_list),n):
            sub_lista = gi_list[i:i+n]
            sub_listas.append(','.join(map(str,sub_lista)))
            dicionario[key] = sub_listas
        
    return dicionario
def buscar_locus(genbank_id):
    #Realiza a procura de novas versões de uma proteína dado um gi obsoleto e 
    #retorna um dicionário contendo o valor de gi antigo relacionado ao valor 
    #de refseq atual.
    #A procura do novo id(Refseq) é dada por meio da existência de um comentário 
    #em cada record retornado no xml buscado na API. Caso uma mesma proteína 
    #tenha sofrido múltiplas atualizações, haverá a busca até se chegar na última versão.
    
    
    
    time.sleep(1)
    handle = Entrez.efetch(db="protein", id=genbank_id, rettype="gb", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    locus_id_novo = {}
    
    for record in records:
        gi_obsoleto = None
        
        if 'GBSeq_comment' in record and record['GBSeq_comment']:
            match_comm = re.search(r'this sequence was replaced by (\w+\.\d+)', record['GBSeq_comment'])
        else:
            match_comm = None
            
        for seqid in record.get("GBSeq_other-seqids", []):
            if seqid.startswith('gi|'):
                gi_obsoleto = seqid
                break

        if match_comm:
            while match_comm:
                refseq_potencial = match_comm.group(1) #.rstrip('.')
                if refseq_potencial.startswith("gi:"):
                    refseq_potencial = refseq_potencial.replace(":", "").lstrip("gi")
                
                time.sleep(1)
                handle = Entrez.efetch(db="protein", id=refseq_potencial, rettype="gb", retmode="xml")
                records = Entrez.read(handle)
                handle.close()

                if 'GBSeq_comment' in records[0] and records[0]['GBSeq_comment']:
                    match_comm = re.search(r'this sequence was replaced by (\w+\.\d+)', records[0]['GBSeq_comment'])
                else:
                    match_comm = None
                    
                info_CDS = []
                if 'GBSeq_feature-table' in records[0] and records[0]['GBSeq_feature-table']:
                    for feature in records[0].get("GBSeq_feature-table",[]):
                        if feature.get("GBFeature_key") == 'CDS':
                            info_CDS.append("FEATURE_LOCATION:"+str(feature.get('GBFeature_location',"NOGBFeature_location")))
                            
                            for qualifier in feature.get('GBFeature_quals',[]):
                                #Qualificador retirado de CDS : 'transl_table': 'NOTransl_table',
                                qualifiers_CDS = {
                                            'gene': 'NOGene',
                                            'locus_tag': 'NOLocus_tag',
                                            'standard_name': 'NOStandard_name',
                                            'coded_by': 'NOCoded_by',
                                            'note': 'NONote',
                                            'db_xref': 'NODb_xref'
                                        }
                                qualifier_name = qualifier.get("GBQualifier_name")
                                
                                if qualifier_name in qualifiers_CDS: 
                                    info_CDS.append(f"{qualifier_name}:"+str(qualifier.get('GBQualifier_value',qualifiers_CDS[qualifier_name])))
                                    
                #atualiza o dicionário com o refseq.version mais recente
                gi_obsoleto = gi_obsoleto or refseq_potencial
                
                locus_id_novo[gi_obsoleto] = ["REFSEQ:"+str(records[0].get("GBSeq_accession-version", refseq_potencial)),
                                              "SEQUENCE:"+str(records[0].get("GBSeq_sequence","")),
                                              "DEFINITION:"+str(records[0].get("GBSeq_definition","NOGBSeq_definition"))
                                              ] + info_CDS
        else:
            #caso não haja substituição, armazena o refseq.version atual
            info_CDS = []
            if 'GBSeq_feature-table' in record and record['GBSeq_feature-table']:
                for feature in record.get("GBSeq_feature-table",[]):
                    if feature.get("GBFeature_key") == 'CDS':
                        info_CDS.append("FEATURE_LOCATION:"+str(feature.get('GBFeature_location',"NOGBFeature_location")))
                        
                        for qualifier in feature.get('GBFeature_quals',[]):
                            #Qualifier CDS retirado: 'transl_table': 'NOTransl_table',
                            qualifiers_CDS = {
                                        'gene': 'NOGene',
                                        'locus_tag': 'NOLocus_tag',
                                        'standard_name': 'NOStandard_name',
                                        'coded_by': 'NOCoded_by',
                                        'note': 'NONote',
                                        'db_xref': 'NODb_xref'
                                    }
                            qualifier_name = qualifier.get("GBQualifier_name")
                            
                            if qualifier_name in qualifiers_CDS: 
                                info_CDS.append(f"{qualifier_name}:"+str(qualifier.get('GBQualifier_value',qualifiers_CDS[qualifier_name])))
                                                       
            gi_obsoleto = gi_obsoleto or str(record.get("GBSeq_accession-version", ""))
            
            locus_id_novo[gi_obsoleto] = ["REFSEQ:"+str(record.get("GBSeq_accession-version", "")),
                                          "SEQUENCE:"+str(record.get("GBSeq_sequence", "NOSEQUENCE")),
                                          "DEFINITION:"+str(record.get("GBSeq_definition","NOGBSeq_definition")),
                                         ] + info_CDS

    return locus_id_novo if records else None
#Menu de escolha e manuseio do programa
while True:
    print("\nMenu de Seleção:")
    print("1 - Cria diretorio (Caso já tenha alguma pasta com o nome\"batches_gi_refseq_atual\" e/ou algum arquivo nela, ele excluirá os mesmos!!!)")
    print("2 - Leitura dos batchs de arquivos e busca dos RefSeq's Atuais.")
    print("0 - Sair")
    
    escolha = input("Por favor, escolha uma opção: \n")
    print("\n")
    
    if escolha == "1":
        #criacao de arquivos .json dos batches
        pasta = "batches_gi_refseq_atual"

        diretorio_atual = os.getcwd()
        caminho_pasta = os.path.join(diretorio_atual, pasta) 

        if os.path.exists(caminho_pasta):
            shutil.rmtree(caminho_pasta)    
        os.mkdir(caminho_pasta)
        print("#####################################################")
        print(f"Criado Diretório:{pasta}")
        print("#####################################################")
        print("\n\n")
    elif escolha == "2":
        inicio = datetime.datetime.now()
        pasta_leitura1 = "kog.txt"
        pasta_leitura2 = "kyva=gb.txt"
        pasta = "batches_gi_refseq_atual"

        diretorio_atual = os.getcwd()
        caminho_pasta = os.path.join(diretorio_atual, pasta)
        caminho_pasta_leitura1 = os.path.join(diretorio_atual, pasta_leitura1)
        caminho_pasta_leitura2 = os.path.join(diretorio_atual, pasta_leitura2)

        #Manuseio dos arquivos kog
        with open (caminho_pasta_leitura1,'r') as file:
            rows = file.readlines()
        kog_dado = [row.strip().split() for row in rows]
        
        #novidade
        kog_descricao = [row for row in kog_dado if len(row)>2]
        pd_Descricao = pd.DataFrame(kog_descricao)

        pd_Descricao['Descrição'] = pd_Descricao.iloc[:,2:19].fillna('').agg(r' '.join, axis=1)
        pd_Descricao = pd_Descricao.drop(columns = pd_Descricao.iloc[:,2:19])
        pd_Descricao = pd_Descricao.rename(columns={0:"OL.FUNC",1:"KOG"})
        #

        kog = pd.DataFrame(kog_dado)
        kog = kog.drop(kog.columns[2:],axis = 1).dropna()
        kog.columns = ['kog_cod_organismos','kog_protein_name']

        with open(caminho_pasta_leitura2,'r') as file:
            rows = file.readlines()
        kyva_gb_dado = [row.strip().split() for row in rows]

        kyva_gb = pd.DataFrame(kyva_gb_dado, columns = ['kog_protein_name','gi'])

        print("#####################################################")   
        print("Leitura dos Arquivos KOG realizada")
        print("#####################################################")
        print("\n\n")


        #Organismos buscados:
        #Cel  Caenorhabditis elegans
        #Dme  Drosophila melanogaster
        #Hsa  Homo sapiens
        #Sce  Saccharomyces cerevisiae
        kog_cod_organismos = ['cel:','dme:','hsa:','sce:']
        
        kog['kog_protein_name'] = kog['kog_protein_name'].apply(limpeza_kog)
        #novidade Criação e atribuição dos KOG
        kog['KOG'] = kog['kog_protein_name'].where(kog['kog_protein_name'].str.match(r'KOG\d+'))
        kog['KOG'] = kog['KOG'].ffill()
        #

        kog_filtrado_proteina_organismos = kog[kog['kog_cod_organismos'].isin(kog_cod_organismos)]
        
        #novidade merge KOG,OL FUNC, DESCRICAO
        kog_filtrado_proteina_organismos = pd.merge(kog_filtrado_proteina_organismos,pd_Descricao,on="KOG",how = 'left')
        #
        
        merge_kog_kyva = pd.merge(kog_filtrado_proteina_organismos, kyva_gb,on='kog_protein_name',how = 'left')
        merge_kog_kyva = merge_kog_kyva.drop_duplicates(subset = ['gi'],keep = 'first')
        
        lista_gi_obsoleto = {}
        for cod in kog_cod_organismos:
            lista_gi_obsoleto[cod] = merge_kog_kyva['gi'][merge_kog_kyva['kog_cod_organismos']==cod].tolist()
        
        #criacao batches para busca com entrez
        lista_batch = cria_batch(lista_gi_obsoleto)
        
        print("#####################################################")
        print("Filtragem dos dados KOG a partir dos organismos modelos realizada")
        print("#####################################################")
        print("\n\n")
        for key in lista_batch.keys():
            print(f"Organismo -> \"{key.replace(':','')}\"")
            
        organismos = input("Digite o código do organismo o qual você quer criar e buscar os refseq's atuais (exatamente como informado acima):")

        organismos = organismos + ':'
        limite_superior = len(lista_batch[organismos])

        #uso do entrez!!!!!
        limite_inferior = int(input(f"Digite um limite inferior para os batches(limite superior: {limite_superior}): "))
        for i in range(limite_inferior,limite_superior):
            inicio_batch = datetime.datetime.now()
            
            buscas_batch_locus = buscar_locus(lista_batch[organismos][i])
            
            #novidade
            info_adicionais = merge_kog_kyva[merge_kog_kyva['kog_cod_organismos']==organismos]

            for gi_key, value in buscas_batch_locus.items():
                
                gi_busca = gi_key.split("|")[1]
                match = info_adicionais[info_adicionais['gi']== gi_busca]
                
                if not match.empty:
                    info_list = match[['kog_protein_name',
                                       'KOG', 'OL.FUNC',
                                       'Descrição']].values.flatten().tolist()
                    for w in range(len(value)):
                        info_list.append(value[w])
                    buscas_batch_locus[gi_key]  = info_list
            #
            
            caminho_pasta = os.path.join(diretorio_atual, pasta)# atual\batches_gi_refseq_atual
            limpeza_cod = organismos.split(":")[0]
            caminho_pasta_organismo = os.path.join(caminho_pasta,f"{limpeza_cod}")#atual\batches_gi_refseq_atual\cel

            if not os.path.exists(caminho_pasta_organismo):
                os.mkdir(caminho_pasta_organismo)
                print("\n")
                print("#####################################################")
                print(f"Criado Diretório: batches_gi_refseq_atual|{limpeza_cod}")
                print("#####################################################")
                print("\n\n")
            
            
            caminho_arquivo = os.path.join(caminho_pasta_organismo,f"batch_{i}.json")
            with open(caminho_arquivo,'w') as arquivo:
                json.dump(buscas_batch_locus,arquivo,indent=4)
            
            #variaveis para plot de tempo
            final_batch = datetime.datetime.now()
            diff_batch = final_batch - inicio_batch
            total_segundos = diff_batch.total_seconds()
            horas = int(total_segundos // 3600)
            minutos = int((total_segundos % 3600) // 60)
            segundos = int(total_segundos % 60)
            fractions = int((total_segundos - int(total_segundos)) * 1000)
            data_hora_formatada = f"{horas:02}:{minutos:02}:{segundos:02}:{fractions:03}"
            print(f"Criado batch_{i}\t||\tTempo Exec. Batch: {data_hora_formatada}\n")
            
        #variaveis para plot de tempo
        
        final = datetime.datetime.now()
        diff = final - inicio
        total_segundos = diff.total_seconds()

        horas = int(total_segundos // 3600)
        minutos = int((total_segundos % 3600) // 60)
        segundos = int(total_segundos % 60)
        fractions = int((total_segundos - int(total_segundos)) * 1000)
        data_hora_formatada = f"{horas:02}:{minutos:02}:{segundos:02}:{fractions:03}"
        print("#####################################################")
        print(f"Tempo total de execução: {data_hora_formatada}")
        print("\n")
        
    elif escolha == "0":
        print("Encerrando o programa...")
        break
    
    else:
        print("Opção inválida. Tente novamente.")
        