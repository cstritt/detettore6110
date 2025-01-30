#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% TSDs to join 5 and 3 prime anchors?
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
   

for side in tsds:
    print('\n')
    for tsd in tsds[side]:
        print(f'{side}\t{tsd}\t{tsds[side][tsd]}')




#%% plot the distribution of cluster sizes
sizes_5 = [len(clusters_5[cl]) for cl in clusters_5]
sizes_3 = [len(clusters_3[cl]) for cl in clusters_3]

# Plot the distribution of cluster sizes
import matplotlib.pyplot as plt
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

