import pandas as pd
import numpy as np
import scipy as sp
from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
import matplotlib.pyplot as plt
from collections import defaultdict
from decimal import Decimal
import math
from k_means_abs_maf import KMeansAbsMaf


class AbsMafAnalyzer:
    def __init__(self, abs_maf_path, detection_power_threshold=0, accession=None, exclude_silent=True, k=None, k_ranges=None):
        self.abs_maf_path = abs_maf_path
        self.detection_power_threshold = detection_power_threshold
        self.exclude_silent = exclude_silent
        self.__read_abs_maf()
        self.k = k
        self.accession = accession
        self.k_ranges = k_ranges if k_ranges else range(1, 15)

    def cluster(self):
        if self.k is None:
            self.__determine_optimal_k()
        else:
            self.__cluster_ccfs(self.k)

    def __determine_optimal_k(self):
        threshold_angle = 7
        k_vs_avg_ssd = {}
        percentage_decrease_at_k = {}
        k_models = {}
        for k in self.k_ranges:
            km = self.__cluster_ccfs(k)
            km_2 = self.__cluster_ccfs(k)
            km_3 = self.__cluster_ccfs(k)

            avg_ssd = (km.ssd['ccf_hat'] + km_2.ssd['ccf_hat'] + km_3.ssd['ccf_hat']) / 3

            k_vs_avg_ssd[k] = avg_ssd
            if k > 1:
                percentage_decrease_at_k[k] = ((k_vs_avg_ssd[k-1] - avg_ssd) / avg_ssd)
            k_models[k] = km

        max_decrease = max(percentage_decrease_at_k.values())
        elbow_k = [k for k, v in percentage_decrease_at_k.items() if v == max_decrease][0]
        best_k = elbow_k
        for k in self.k_ranges[elbow_k:]:
            radian_angle = np.arctan(k_vs_avg_ssd[k] - k_vs_avg_ssd[k-1])
            degrees_angle = 180 * radian_angle / math.pi
            print(k, degrees_angle)
            if -degrees_angle < threshold_angle:
                best_k = k
                break

        plt.plot(k_vs_avg_ssd.keys(), k_vs_avg_ssd.values())
        plt.xlabel('Number of Clusters (k)')
        plt.ylabel('Sum of Squares Distance')
        plt.title('Sum of Squares Distance vs. Number of Clusters')

        plt.axvline(best_k, color='g', linestyle='dashed', linewidth=2)

        self.__plot_final_clusters(k_models[best_k])

    def __cluster_ccfs(self, num_clusters):
        # Cluster
        km = KMeansAbsMaf(self.abs_maf, columns=['ccf_hat'])
        km.cluster(num_clusters)
        return km

    def __plot_final_clusters(self, km):
        num_clusters = len(km.centers)
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

        plt.title('Distribution of Cancer Cell Fraction for SNPs {}'
                  .format('in {}'.format(self.accession) if self.accession else None))
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