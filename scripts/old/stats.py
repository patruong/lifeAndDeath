#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:34:38 2020

@author: ptruong
"""

import numpy as np
import pandas as pd

def count_nan(df):
    count_nan = len(df) - df.count()
    return count_nan


def normalize(df):
    """
    mean normalization column-wise
    """
    normalized_df=(df-df.mean())/df.std()
    return normalized_df