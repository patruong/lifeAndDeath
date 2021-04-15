#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:48:58 2021

@author: ptruong
"""

def drop_zero_row(df):
    return df[(df.T != 0).any()]