#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:34:32 2018

@author: qinggang
"""
import numpy as np
import matplotlib.pylab as plt
from nilearn import plotting

def plot_matrices(matrices, matrix_kind):
    n_matrices = len(matrices)
    fig = plt.figure(figsize = (n_matrices * 4, 4))
    for n_sub, matrix in enumerate(matrices):
        plt.subplot(1, n_matrices, n_sub + 1)
        matrix = matrix.copy()
        np.fill_diagonal(matrix, 0)
        vmax = np.max(np.abs(matrix))
        title = '{0}, subject {1}'.format(matrix_kind, n_sub)
        plotting.plot_matrix(matrix, vmin = -vmax, vmax = vmax,
                             cmap = 'RdBu_r', title = title, figure = fig,
                             colorbar = False)
        
        
plot_matrices(corr_matrix[:4], 'correlation')