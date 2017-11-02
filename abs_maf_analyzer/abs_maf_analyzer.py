import pandas as pd
import numpy as np
import scipy as sp
from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
import matplotlib.pyplot as plt
from collections import defaultdict
from decimal import Decimal
from k_means_abs_maf import KMeansAbsMaf


class AbsMafAnalyzer:
    def __init__(self, abs_maf_path, detection_power_threshold=0, exclude_silent=True):
        self.abs_maf_path = abs_maf_path
        self.detection_power_threshold = detection_power_threshold
        self.exclude_silent = exclude_silent
        self.__read_abs_maf()

    def cluster_ccfs(self):
        # Cluster
        num_clusters = 4
        km = KMeansAbsMaf(self.abs_maf, columns=['ccf_hat'])
        km.cluster(4)

        groups = []
        for i in range(num_clusters):
            groups.append([])

        for index, row in self.abs_maf.iterrows():
            group_index = km.clusters[index]
            groups[group_index].append({'ccf_hat': row.ccf_hat, 'dp': row.detection_power})

        # Create plot
        fig = plt.figure()
        fig.set_size_inches(10, 4)
        ax = fig.add_subplot(1, 1, 1, facecolor="1.0")

        centroid_values = [vector[0] for vector in km.centers.values]
        for data, color, label, centroid in zip(groups,
                                                ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'),
                                                range(num_clusters),
                                                centroid_values):
            cluster_size = len(data)
            # Change size of centroid depiction depending on number of points at value
            unit_size = 20
            fraction_in_cluster = cluster_size/len(self.abs_maf)
            centroid_marker_size = fraction_in_cluster*unit_size + 8

            ax.set_ylim(bottom=-.2, top=1.1)
            ax.set_yticks(np.linspace(0, 1, 11))

            ax.set_xlabel("Cancer Cell Fraction")
            ax.set_ylabel("Detection Power")
            ax.plot([centroid], [-.1], '.', c='k', markeredgewidth=0, markerfacecolor=color,
                    markeredgecolor='k', markersize=centroid_marker_size)

            ax.scatter([d['ccf_hat'] for d in data], [d['dp'] for d in data],
                       alpha=1, c=color, edgecolors='none', s=unit_size, label=label)

        plt.title('Distribution of Cancer Cell Fraction for SNPs')
        plt.axhspan(-.2, 0, facecolor='0.2', alpha=0.4)
        plt.show()

    def __read_abs_maf(self):
        abs_maf = pd.read_csv('{}'.format(self.abs_maf_path), sep='\t')
        abs_maf = abs_maf.loc[:, ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position',
                                  'Variant_Classification', 'Variant_Type', 'Reference Alelle',
                                  'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'ccf_hat', 'detection_power']]

        if self.exclude_silent:
            abs_maf = abs_maf[abs_maf.Variant_Classification != 'Silent']
        self.abs_maf = abs_maf[abs_maf.detection_power >= self.detection_power_threshold]