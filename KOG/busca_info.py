# -*- coding: utf-8 -*-

import datetime
import pandas as pd
from Bio import Entrez
import os
import json
import shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import time


#funcao
def buscar_info(batch_dict):
    dict_refseq_sequencia = {}
    for proteina in batch_dict.values():
        dict_refseq_sequencia[proteina[4]] = [proteina[6],proteina[5]]
    
    return dict_refseq_sequencia

while True:
    print("\nMenu de Seleção:")
    print("1 - Cria diretorio (Caso já tenha alguma pasta com o nome \"batches_sequencia\" e/ou algum arquivo nela, ele excluirá os mesmos!!!)")
    print("2 - Leitura dos batchs de arquivos batch_X.json da Pasta \"batches_gi_refseq_atual\"")
    print("3 - cria o arquivo .fasta de cada batch")
    print("4 - Junção dos arquivos fasta e adição e criação arquivo .csv")
    print("0 - Sair")
    
    escolha = input("Por favor, escolha uma opção:\n")
    if escolha == "1":
        #criação diretorio
        pasta = "batches_sequencias"
        diretorio_atual = os.getcwd()
        caminho_pasta = os.path.join(diretorio_atual, pasta) 

        if os.path.exists(caminho_pasta):
            shutil.rmtree(caminho_pasta)    
        os.mkdir(caminho_pasta)
        print("#####################################################")
        print("Criado Diretório")
        print("#####################################################")
        print('\n')
    elif escolha =="2":
        
        
        pasta_leitura = "batches_gi_refseq_atual"
        
        lista_de_dicionarios = []
        diretorio_atual = os.getcwd()
        caminho_pasta_leitura = os.path.join(diretorio_atual, pasta_leitura)
        arquivos_json = [arquivo for arquivo in os.listdir(caminho_pasta_leitura)]
        
        for pasta in arquivos_json:
            print(f"Pasta -> \"{pasta}\"")
        
        organismo = input("Digite a pasta do organismo o qual você quer buscar as sequencias (exatamente como informado acima):\n")
        
        caminho_pasta_organismos = os.path.join(caminho_pasta_leitura,organismo)
        arquivos_json = [arquivo for arquivo in os.listdir(caminho_pasta_organismos) if arquivo.endswith(".json")]
        arquivos_json.sort(key=lambda x: int(x.split('_')[1].split('.')[0]))
        
        limite_superior = len(arquivos_json)
        m = int(input(f"Digite um limite inferior para os batches(limite superior: {limite_superior}): // 0 - todos arquivos: \n"))#variavel para inicio da iteração dos batchs
        print('\n')
        #batch_ipg_200
        for i in range(m,len(arquivos_json)):
             # Caminho completo do arquivo
             nome_do_arquivo = arquivos_json[i]
             caminho_do_arquivo = os.path.join(caminho_pasta_organismos, nome_do_arquivo)
             # Abre o arquivo e carrega o conteúdo como um dicionário
             with open(caminho_do_arquivo, 'r') as arquivo:
                 dicionario = json.load(arquivo)
                 lista_de_dicionarios.append(dicionario)     
             print(f"Arquivo {nome_do_arquivo} foi lido e adicionado à lista.")
        print('\n')      
    elif escolha =="3":
        inicio = datetime.datetime.now()
        
        pasta_escrita = "batches_sequencias"
        diretorio_atual = os.getcwd()
        caminho_pasta_escrita = os.path.join(diretorio_atual, pasta_escrita)
        try:
            if not organismo:
                raise ValueError("Variavel para Organismo nao foi definida no passo 2, por favor execute o passo 2 e depois o 3 passo.\n")
                caminho_pasta_organimo_escrita = os.path.join(caminho_pasta_escrita, organismo)
        except Exception as e:
            print('\n')
            print(f"Ocorreu um erro: {e}")
            print('\n')
            print("Variável para Organismo não foi definida no passo 2, por favor execute o passo dois e depois passo três.\n")
            print('\n')
            print("Encerrando o programa...")
            sys.exit(1)
        caminho_pasta_organismo_escrita = os.path.join(caminho_pasta_escrita, organismo)
        print("#####################################################")
        if not os.path.exists(caminho_pasta_organismo_escrita):
            os.mkdir(caminho_pasta_organismo_escrita)
            print("\n")
            print("#####################################################")
            print(f"Criado Diretório: batches_sequencias|{organismo}")
            print("#####################################################")
            print("\n\n")
        #chgamada escolha 3
        limite_superior = len(lista_de_dicionarios)
        limite_inferior = int(input(f"Digite um limite inferior para os batches(limite superior: {limite_superior}):\n"))
        for i in range(limite_inferior,len(lista_de_dicionarios)):
            inicio_batch = datetime.datetime.now()
            dicionario_sequencias = buscar_info(lista_de_dicionarios[i])
            
            caminho_arquivo_escrita = os.path.join(caminho_pasta_organismo_escrita,f"batch_sequencia_{i}.fasta")
            
            records = []
            
            for refseq, (descricao,sequencia) in dicionario_sequencias.items():
                record = SeqRecord(Seq(sequencia.split(':',1)[1]), id=refseq.split(':',1)[1], description=descricao.split(':',1)[1])
                records.append(record)
            with open(caminho_arquivo_escrita,'w') as arquivo:
                #formato fasta
                SeqIO.write(records,arquivo,"fasta")
            final_batch = datetime.datetime.now()
            diff_batch = final_batch - inicio_batch
            total_segundos_ = diff_batch.total_seconds()
            horas_ = int(total_segundos_ // 3600)
            minutos_ = int((total_segundos_ % 3600) // 60)
            segundos_ = int(total_segundos_ % 60)
            fractions_ = int((total_segundos_ - int(total_segundos_)) * 1000)
            data_hora_formatada_ = f"{horas_:02}:{minutos_:02}:{segundos_:02}:{fractions_:03}"
            print(f"Foi Criado o batch_sequencia_{i}.fasta\t||\tTempo Exec. Batch: {data_hora_formatada_}\n")
        #variaveis para plot de tempo
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
    elif escolha =="4":
        pasta_sequencias = "batches_sequencias"
        
        diretorio_atual = os.getcwd()
        caminho_pasta_leitura_sequencia = os.path.join(diretorio_atual, pasta_sequencias)
        pastas = [arquivo for arquivo in os.listdir(caminho_pasta_leitura_sequencia)]
        
        for pasta in pastas:
            caminho_leitura_fasta =  os.path.join(caminho_pasta_leitura_sequencia, pasta)
            arquivos_fasta = [arquivo for arquivo in os.listdir(caminho_leitura_fasta) if arquivo.endswith(".fasta")]
            arquivos_fasta.sort(key=lambda x: int(x.split('_')[2].split('.')[0]))
        
            pasta_final = 'resultado'
        
            join_pastafinal = os.path.join(diretorio_atual,pasta_final)
            if not os.path.exists(join_pastafinal):
                os.mkdir(join_pastafinal)
                print("\n")
                print("#####################################################")
                print(f"Criado Diretório: {pasta_final}")
                print("#####################################################")
                print("\n\n")
            nome_arquivo = os.path.join(join_pastafinal, f'{pasta_final}_{pasta}.fasta')
            registros_organismo = []
            for arquivo in arquivos_fasta:
                caminho_arquivo = os.path.join(caminho_leitura_fasta, arquivo)
                with open(caminho_arquivo, "r") as input_handle:
                    registros = list(SeqIO.parse(input_handle, "fasta"))
                    registros_organismo.extend(registros)
            with open(nome_arquivo, "w") as output_handle:
                SeqIO.write(registros_organismo, output_handle, "fasta")
                            
            print(f"Arquivo final criado: {pasta_final}_{pasta}.fasta\n")
        #CSV
        inicio = datetime.datetime.now()
        pasta_json = "batches_gi_refseq_atual"
        diretorio_atual = os.getcwd()
        caminho_pasta_leitura_info = os.path.join(diretorio_atual,pasta_json)
        print("\n")
        print("#####################################################")
        print("Criando arquivo csv. Espere...")
        print("#####################################################")
        print("\n\n")
        pastas =  [arquivo for arquivo in os.listdir(caminho_pasta_leitura_info)]
        lista_info = []

        for pasta in pastas:  
            caminho_leitura_info =  os.path.join(caminho_pasta_leitura_info, pasta)
            arquivos_json = [arquivo for arquivo in os.listdir(caminho_leitura_info) if arquivo.endswith(".json")]
            arquivos_json.sort(key=lambda x: int(x.split('_')[1].split('.')[0]))
            for arquivo in arquivos_json:
                caminho_arquivo = os.path.join(caminho_leitura_info, arquivo)
                with open(caminho_arquivo, "r") as info_json:
                    lista_informacao = json.load(info_json)
                    lista_info.append(lista_informacao) 

        df_csv = pd.DataFrame(columns=[
                                'KOG_GI',
                                'KOG_PROTEIN_NAME',
                                'KOG', 
                                'KOG_OL.FUNC',
                                'KOG_DESCRICAO',
                                'NCBI_REFSEQ_ACC_VERSION',
                                'NCBI_SEQUENCE',
                                'NCBI_DEFINITION',
                                'NCBI_FEATURE_LOCATION',
                                'NCBI_GENE',
                                'NCBI_LOCUS_TAG',
                                'NCBI_STANDARD_NAME',
                                'NCBI_CODED_BY',
                                'NCBI_NOTE',
                                'NCBI_DB_XREF'
                                ]
                            )
        qualifiers_CDS = {
                    'REFSEQ':'NCBI_REFSEQ_ACC_VERSION',
                    'SEQUENCE':'NCBI_SEQUENCE',
                    'DEFINITION':'NCBI_DEFINITION',
                    'FEATURE_LOCATION':'NCBI_FEATURE_LOCATION',
                    'gene':'NCBI_GENE',
                    'locus_tag': 'NCBI_LOCUS_TAG',
                    'standard_name': 'NCBI_STANDARD_NAME',
                    'coded_by': 'NCBI_CODED_BY',
                    'note': 'NCBI_NOTE',
                    'db_xref': 'NCBI_DB_XREF'
                }
           
        for entrada in lista_info:
            for k,v in entrada.items():
                data = {
                           'KOG_GI': k,
                           'KOG_PROTEIN_NAME':v[0],
                           'KOG':v[1], 
                           'KOG_OL.FUNC':v[2],
                           'KOG_DESCRICAO':v[3],
                           }
                for l in range(4,len(v)):
                    split = v[l].split(':',1)
                    feature = split[0]
                    valor = split[1] if len(split) > 1 else v[l]
                    if feature in qualifiers_CDS:
                        chave = qualifiers_CDS[feature]
                        if chave in data:
                            if not isinstance(data[chave],list):
                                data[chave] = [data[chave]]
                            data[chave].append(valor)
                        else:
                            data[chave] = valor
                df_csv.loc[len(df_csv.index)] = data
                           
        #ADICAO ID UNIPROT
        txt_uniprot = 'gene_refseq_uniprotkb_collab.txt'
        chunksize = 10**6 #Alterar o valor caso esteja com dificuldades em processar o dado (chunks menores utilizaram menos memória porém processo dura mais)
        df_csv['UNIPROTKB_PROTEIN_ACC'] = None

        NCBI_REFSEQ_ACCVER = set(df_csv['NCBI_REFSEQ_ACC_VERSION'])

        path_uniprot = os.path.join(diretorio_atual,txt_uniprot)

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
                df_csv.loc[df_csv['NCBI_REFSEQ_ACC_VERSION'] == ncbi_acc,'UNIPROTKB_PROTEIN_ACC'] = uniprot_acc
        
        nome_csv= os.path.join(join_pastafinal, 'resultado_csv.csv')
        df_csv.to_csv(nome_csv,header=True,index=False)
        
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