#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 03:13:49 2020

@author: ptruong
"""

import requests

import pandas as pd

def go_class(code):
    code = code
    data = requests.get("https://www.uniprot.org/uniprot/"+ code + ".txt").content.decode("utf-8").split("\n")
    
    # Gene-Ontology protein Classification
    data_ = []
    for i in data:
        vals = i.split()
        try:
            if vals[1] == "GO;":
                print(vals)
        except:
            continue
        #data_.append(i.split())

go_class("A0A0D9SFH9")
go_class("A0A075B6E5")
go_class("P13726")
go_class("P13726-2")
go_class("Q99986")
go_class("H0YJ50")


    

df = pd.DataFrame(data)

df_logFc.columns[10000] # from lifeAndDeath.py



