# -*- coding: utf-8 -*-
"""
Recebe o arquivo de  essential_genes e separa por organismos assim como
retira as entradas onde não se tem sequência

Esse código pode criar o arquivo de sequência das proteínas ou dos genes

Ele cria arquivos no formato fasta para cada organismo da coluna Organism
"""

import pandas as pd
import os

dict_final = {}
visualizer = pd.DataFrame()

diretorio_atual = os.getcwd()
df = pd.read_csv(os.path.join(diretorio_atual, "essential_genes.csv"), sep = ",")
 
    
for organism, group in df.groupby("Organism"):
    fasta_path = os.path.join(diretorio_atual, f"gene_{organism}.fasta")
    
    with open(fasta_path, "w") as fasta_file:
        for _, row in group.iterrows():
            funcao = row["Function"].strip('"')
            id_proteina = row["RefSeq"]
            #Mudar para seq_Prot se quiser a sequência da proteína
            sequencia = row["Seq_Gene"]
            # Escreve no formato de arquivo FASTA
            if sequencia != 'Notavailablenow.':
                fasta_file.write(f">{id_proteina} {funcao} [{organism}]\n{sequencia}\n")
    
    print(f"FASTA file created for organism: {organism} at {fasta_path}")