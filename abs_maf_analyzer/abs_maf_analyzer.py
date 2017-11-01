import pandas as pd
import numpy as np
import scipy as sp
from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
import matplotlib.pyplot as plt
from k_means_plus_plus import *

class AbsMafAnalyzer:
    def __init__(self, abs_maf_path):
        self.abs_maf_path = abs_maf_path
        self.__read_abs_maf()

    def cluster_ccfs(self):
        # Cluster
        kmpp = KMeansPlusPlus(self.abs_maf, 3, columns=['ccf_hat'])
        kmpp.cluster()
        print(kmpp.clusters)
        # Get a scatterplot that's color-coded by cluster
        #colors = [
        #    "red" if x == 0 else "blue" if x == 1 else "green" for x in kmpp.clusters]
        #plt.scatter(data['x'], data['y'], s=5, c=colors)
        #plt.savefig("three_clusters_clusters.png")

    def __read_abs_maf(self):
        abs_maf = pd.read_csv('{}'.format(self.abs_maf_path), sep='\t')
        abs_maf = abs_maf.loc[:, ['Hugo_Symbol', 'Chromosome', 'Start_position', 'End_position',
                                  'Variant_Classification', 'Variant_Type', 'Reference Alelle',
                                  'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'ccf_hat', 'detection_power']]
        self.abs_maf = abs_maf[abs_maf.Variant_Classification != 'Silent']
        pass