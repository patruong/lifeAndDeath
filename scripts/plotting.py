#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:45:55 2020

@author: ptruong
"""

def create_markers(cell_lines, states):
    # Create marker dict
    marker_dict = dict()
    fillstyles = []
    markers = ["o","1",
               "^","x",
               "s","+"]
    count = 0
    for i in np.unique(cell_lines):
        for j in np.unique(states):
            marker_for = i+"_"+j        
            marker_dict.update({marker_for: markers[count]})
            count += 1
    return market_dict
    

