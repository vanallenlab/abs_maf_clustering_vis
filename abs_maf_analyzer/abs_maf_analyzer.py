import pandas as pd
import numpy as np
import scipy as sp
from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
import matplotlib.pyplot as plt
from collections import defaultdict
from k_means_plus_plus import *


class AbsMafAnalyzer:
    def __init__(self, abs_maf_path):
        self.abs_maf_path = abs_maf_path
        self.__read_abs_maf()

    def cluster_ccfs(self):
        # Cluster
        num_clusters = 5
        kmpp = KMeansPlusPlus(self.abs_maf, num_clusters, columns=['ccf_hat'])
        kmpp.cluster()
        print(kmpp.clusters)
        centroid_values = [vector[0] for vector in kmpp.centers.values]
        groups = []
        for i in range(num_clusters):
            groups.append([])

        for index, row in self.abs_maf.iterrows():
            group_index = kmpp.clusters[index]
            groups[group_index].append(row.ccf_hat)

        colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w')
        group_labels = range(num_clusters)

        # Create plot
        fig = plt.figure()
        fig.set_size_inches(10, 4)
        ax = fig.add_subplot(1, 1, 1, facecolor="1.0")

        for data, color, label, centroid in zip(groups, colors, group_labels, centroid_values):
            # Change size depending on number of points at value
            unit_size = 40

            # Remove duplicates
            data_deduped = list(set(data))

            cluster_size = len(data)

            ax.set_ylim(bottom=-.0002, top=.0005)
            ax.set_yticklabels([])
            ax.set_xlabel("Cancer Cell Fraction")
            ax.plot([centroid], [.0003], '.', c='k', markeredgewidth=1, markerfacecolor=color, alpha=0.5,
                    markeredgecolor='k', markersize=(cluster_size/len(self.abs_maf)*unit_size)+5)
            ax.annotate("n = {}".format(len(data)), (centroid, .0003))
            ax.scatter(data_deduped, np.zeros(len(data_deduped)), alpha=1, c=color, edgecolors='none', s=unit_size, label=label)

        plt.title('Distribution of Cancer Cell Fraction for SNPs')
        plt.show()

    def __read_abs_maf(self):
        abs_maf = pd.read_csv('{}'.format(self.abs_maf_path), sep='\t')
        abs_maf = abs_maf.loc[:, ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position',
                                  'Variant_Classification', 'Variant_Type', 'Reference Alelle',
                                  'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'ccf_hat', 'detection_power']]
        self.abs_maf = abs_maf[abs_maf.Variant_Classification != 'Silent']
        pass