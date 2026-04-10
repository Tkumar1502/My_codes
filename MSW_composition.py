# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 17:18:51 2025

@author: yoelr
"""
import numpy as np
import pandas as pd
import biosteam as bst
import seaborn as sns
import matplotlib.pyplot as plt
from plastics import strap
import os

folder = os.path.dirname(__file__)
folder = os.path.join(folder, 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')

def plot_MRF_residue_composition():
    fs = 12
    pm = strap.STRAPMSWProcess(preprocessing=True)
    MRF_residue = pm.feedstock
    IDs = [
        'BiogenicMaterial', 
        'PEPP',
        'BulkPlastic', 
        'NonResidue',
    ]
    values = MRF_residue.imass[IDs]
    values /= values.sum()
    values *= 100
    IDs[-1] = 'Metal/Glass/Rock'
    df = pd.DataFrame(values, index=IDs, columns=['MRF residues'])
    df = df.sort_values(by=['MRF residues'], ascending=False)
    sns.set(style='ticks')
    bst.set_figure_size(aspect_ratio=1.6, width=6.6142 / 3.5)
    bst.set_font(size=fs)
    fig, ax = plt.subplots(1, 1)
    plt.sca(ax)
    df.T.plot.bar(stacked=True, ax=ax, fontsize=fs)
    plt.ylabel(r'Composition [% wt]', fontsize=fs)
    ax = plt.gca()
    plt.yticks([0, 25, 50, 75, 100])
    ax.tick_params(axis='x', which='major', length=6,
                   direction="inout", rotation=0)
    ax.tick_params(axis='y', which='major', length=6,
                   direction="inout")
    ax.get_legend().remove()
    ax.spines[['right', 'top']].set_visible(False)
    plt.subplots_adjust(left=0.35)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MRF_residue_composition.{i}')
        plt.savefig(file, dpi=900, transparent=True)
    
if __name__ == '__main__':
    plot_MRF_residue_composition()