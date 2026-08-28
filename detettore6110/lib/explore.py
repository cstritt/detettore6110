#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Exploratory analysis functions for detettore6110 results.
This module contains functions for exploring and visualizing results.
"""

def analyze_tsds(cluster_d):
    """
    Analyze TSDs to join 5 and 3 prime anchors.
    
    Parameters
    ----------
    cluster_d : dict
        Dictionary containing cluster data for both sides
        
    Returns
    -------
    dict
        Dictionary with TSD counts per side
    """
    tsds = {}
    for side in cluster_d:
        tsds[side] = {}
        
        for cluster_id in cluster_d[side]:
            if side == '5':
                tsd = cluster_d[side][cluster_id].consensus.seq[-3:]
            elif side == '3':
                tsd = cluster_d[side][cluster_id].consensus.seq[:3]
            
            if tsd not in tsds[side]:
                tsds[side][tsd] = 0
            tsds[side][tsd] += 1
    
    return tsds


def plot_cluster_sizes(clusters_5, clusters_3):
    """
    Plot the distribution of cluster sizes.
    
    Parameters
    ----------
    clusters_5 : dict
        Dictionary of 5' prime clusters
    clusters_3 : dict
        Dictionary of 3' prime clusters
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed. Skipping plots.")
        return
    
    sizes_5 = [len(clusters_5[cl]) for cl in clusters_5]
    sizes_3 = [len(clusters_3[cl]) for cl in clusters_3]

    # Plot the distribution of cluster sizes
    plt.hist(sizes_5)
    plt.xlabel('Cluster size')
    plt.ylabel('Number of clusters')
    plt.title('Cluster sizes for 5\' side')
    plt.show()
    plt.hist(sizes_3)
    plt.xlabel('Cluster size')
    plt.ylabel('Number of clusters')
    plt.title('Cluster sizes for 3\' side')
    plt.show()
