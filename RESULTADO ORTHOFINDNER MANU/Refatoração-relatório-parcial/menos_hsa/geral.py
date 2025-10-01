# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 23:15:49 2025

@author: Emanu
"""

import pandas as pd
import os

proteinas_man = {}
visualizer = pd.DataFrame()

diretorio_atual = os.getcwd()


for arquivo_merge in os.listdir(diretorio_atual):
    nome = arquivo_merge
    
    
