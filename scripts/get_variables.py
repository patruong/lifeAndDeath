#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:43:28 2021

@author: ptruong
"""

def get_cell_line_states_replicates_from_reporter_intensity_cols(reporter_intensity_corrected_cols):
    cell_lines, states, replicates = [], [], []
    for x in reporter_intensity_corrected_cols:
        cell_line = x.split(" ")[-1].split("_")[0]
        state = x.split(" ")[-1].split("_")[1]
        replicate = x.split(" ")[-1].split("_")[2]
        cell_lines.append(cell_line), states.append(state), replicates.append(replicate)
    cell_lines = list(dict.fromkeys(cell_lines).keys())
    states = list(dict.fromkeys(states).keys())
    replicates = list(dict.fromkeys(replicates).keys())
    return cell_lines, states, replicates    


def get_cell_lines_states_replicates():
    cell_lines = ['A549', 'MCF7', 'RKO']
    states = ['D', 'S']
    replicates = ['Rep1', 'Rep2', 'Rep3']
    return cell_lines, states, replicates