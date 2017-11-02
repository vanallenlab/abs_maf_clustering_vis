from pandas import DataFrame, Series
import numpy as np
from numbers import Integral
import copy


class KMeansAbsMaf:
    def __init__(self, data_frame, columns=None, max_iterations=None):
        if not isinstance(data_frame, DataFrame):
            raise Exception("data_frame argument is not a pandas DataFrame")
        elif data_frame.empty:
            raise Exception("The given data frame is empty")
        if max_iterations is not None and max_iterations <= 0:
            raise Exception("max_iterations must be positive!")

        self.data_frame = data_frame  # m x n
        self.numRows = data_frame.shape[0]  # m

        # k x n, the i,j entry being the jth coordinate of center i
        self.centers = None

        # m x k , the i,j entry represents the distance
        # from point i to center j
        # (where i and j start at 0)
        self.distance_matrix = None
        self.clusters = None
        self.max_iterations = max_iterations
        self.k = 0

        if columns is None:
            self.columns = data_frame.columns
        else:
            for col in columns:
                if col not in data_frame.columns:
                    raise Exception(
                        "Column '%s' not found in the given DataFrame" % col)
                if not self.__is_numeric(col):
                    raise Exception(
                        "The column '%s' is either not numeric or contains NaN values" % col)
            self.columns = columns

    def _populate_initial_centers(self):
        """Randomly choose k initial centroids, ensuring that they all have a different value"""
        rows = list()

        while len(rows) < self.k:
            dice_roll = np.random.rand()
            index = int(round(dice_roll * self.numRows, 0)-1)
            point = self.data_frame[self.columns].iloc[index, :]
            duplicate_point = False
            for existing_point in rows:
                if point[0] == existing_point[0]:
                    duplicate_point = True
            if not duplicate_point:
                point_copy = copy.copy(point)
                point_copy = point_copy.rename(len(rows))
                rows.append(point_copy)

        self.centers = DataFrame(rows, columns=self.columns)

    def __compute_distances(self):
        """Compute distances from each of k centroids for all points"""
        column_dict = {}
        for centroid_index in range(self.k):
            column_dict[centroid_index] = self.__distances_from_point(self.centers.iloc[centroid_index, :])
        self.distance_matrix = DataFrame(column_dict, columns=range(self.k))

    def __get_clusters(self):
        """Compute closest centroid for each point"""
        index_of_closest_centroid = [np.argmin(self.distance_matrix.iloc[i, :], axis=1) for i in range(self.numRows)]
        self.clusters = Series(index_of_closest_centroid, index=self.data_frame.index)

    def __compute_new_centers(self):
        if self.centers is None:
            raise Exception("Centers not initialized!")

        if self.clusters is None:
            raise Exception("Clusters not computed!")

        for i in range(self.k):
            cluster_mean = self.data_frame[self.columns].ix[self.clusters == i].mean()
            self.centers.set_value(i, self.columns, cluster_mean)

    def __compute_ssd(self):
        """Compute sum of square distance between cluster centroid and members over all clusters"""
        total_ssd = 0
        for i in range(self.k):
            centroid_value = - self.centers.iloc[i]
            ssd = np.power(self.data_frame[self.columns].ix[self.clusters == i] - centroid_value, 2).sum()
            total_ssd += ssd
        self.ssd = total_ssd

    def cluster(self, k):
        if not isinstance(k, Integral) or k <= 0:
            raise Exception("k must be a positive integer")
        else:
            self.k = k

        self._populate_initial_centers()
        self.__compute_distances()
        self.__get_clusters()
        self.__compute_ssd()

        counter = 0

        while True:
            counter += 1
            previous_clusters = self.clusters.copy()

            self.__compute_new_centers()
            self.__compute_distances()
            self.__get_clusters()
            self.__compute_ssd()

            if self.max_iterations is not None and counter >= self.max_iterations:
                break
            elif all(self.clusters == previous_clusters):
                break

    def __distances_from_point(self, point):
        """Calculate sum of square distances between given point and all other points"""
        return np.power(self.data_frame[self.columns] - point, 2).sum(axis=1)

    def __is_numeric(self, col):
        return all(np.isreal(self.data_frame[col])) and not any(np.isnan(self.data_frame[col]))
