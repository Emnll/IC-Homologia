# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 18:12:08 2025

@author: Jean
"""

import os
import re
import pandas as pd

pasta_base = os.getcwd()

arquivo_cog = os.path.join(pasta_base,"Resultado_cog_GO_left.csv")
arquivo_kog = os.path.join(pasta_base,"Resultado_kog_GO_left.csv")
arquivo_costa = os.path.join(pasta_base,"Resultado_costa_GO_left.csv") 

#putative attribute -> GO_PROTEIN_NAME(col:30)
df_cog = pd.read_csv(arquivo_cog,
            sep=",", low_memory=False)
#putative attribute -> NCBI_DEFINITION(col:8) GO_PROTEIN_NAME(col:23)
df_kog = pd.read_csv(arquivo_kog,
            sep=",")
#putative attribute -> function(col:10), GO_protein_name(col:22)
df_costa = pd.read_csv(arquivo_costa,
            sep=",")    


def putative(x):
    p = re.search(r".*[Pp]utativ\w*.*",str(x)) is not None
    if(p):
        return True
    else:
        return False
def putative_mask(df):
    #pode ser feito assim porém com menos desempenho em casos de dataframe maior
    #df_kog = df_kog[~df_kog.applymap(lambda x: bool(re.search(r"[Pp]utativ\w*", str(x))))]
    
    mask = df.astype(str).apply(lambda col: col.str.contains(r"[Pp]utativ\w*", regex=True, na=False))
    return mask

#mask = df.astype(str).apply(lambda col: col.str.contains(r"[Pp]utativ\w*", regex=True, na=False))
#df_clean = df[~mask.any(axis=1)]

#cog 40589 - 40584 -> 5 entradas putativas
#costa 4285 - 4233 -> 52 entradas putativas
#kog 22470 - 21552 -> 918 entradas putativas

df_cog = df_cog[~putative_mask(df_cog).any(axis=1)]
df_kog = df_kog[~putative_mask(df_kog).any(axis=1)]
df_costa = df_costa[~putative_mask(df_costa).any(axis=1)]

nome_arquivo_cog = os.path.join(pasta_base,"Resultado_cog_GO_left_sem_putative.csv")
nome_arquivo_kog = os.path.join(pasta_base,"Resultado_kog_GO_left_sem_putative.csv")
nome_arquivo_costa = os.path.join(pasta_base,"Resultado_costa_GO_left_sem_putative.csv")

df_cog.to_csv(nome_arquivo_cog,index=False)
df_kog.to_csv(nome_arquivo_kog,index=False)
df_costa.to_csv(nome_arquivo_costa,index=False)