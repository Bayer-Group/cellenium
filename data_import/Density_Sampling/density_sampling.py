# !/usr/bin/env python


# Density_Sampling/Density_Sampling.py

# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com, ggiecold@jimmy.harvard.edu

# slight modifications (mostly for Python 3) by Carsten Jahn


"""For a data-set comprising a mixture of rare and common populations,
density sampling gives equal weights to selected representatives
of those distinct populations.

Density sampling is a balancing act between signal and noise. Indeed, while
it increases the prevalence of rare populations, it also increases the prevalence
of noisy sample points that would happen to have their local densities larger than
an outlier density computed by Density_Sampling.

An illustration of how to use the present module is in order:

>>> iris = datasets.load_iris()
>>> Y = iris.target
>>> X_reduced = PCA(n_components = 3).fit_transform(iris.data)

>>> plot_PCA(X_reduced, Y, 'the whole Iris data-set')

>>> sampled_indices = density_sampling(X_reduced, metric = 'euclidean', desired_samples = 50)
>>> downsampled_X_reduced = X_reduced[sampled_indices, :]
>>> downsampled_Y = Y[sampled_indices]

>>> plot_PCA(downsampled_X_reduced, downsampled_Y, 'the Iris data-set\ndown-sampled to about 50 samples')

Reference
---------
Giecold, G., Marco, E., Trippa, L. and Yuan, G.-C.,
"Robust Lineage Reconstruction from High-Dimensional Single-Cell Data".
ArXiv preprint [q-bio.QM, stat.AP, stat.CO, stat.ML]: http://arxiv.org/abs/1601.02748
"""

import numbers
import numpy as np
import operator
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import radius_neighbors_graph
from sys import exit
from tempfile import NamedTemporaryFile
from functools import reduce

__all__ = ['get_local_densities', 'density_sampling']


def get_chunk_size(N, n):
    mem_free = 10000000  # just assume 10GB here
    chunk_size = int(((mem_free - 1400000) * 1000) / (4 * n * N))
    return chunk_size


def median_min_distance(data, metric):
    """This function computes a graph of nearest-neighbors for each sample point in
        'data' and returns the median of the distribution of distances between those
        nearest-neighbors, the distance metric being specified by 'metric'.

    Parameters
    ----------
    data : array of shape (n_samples, n_features)
        The data-set, a fraction of whose sample points will be extracted
        by density sampling.

    metric : string
        The distance metric used to determine the nearest-neighbor to each data-point.
        The DistanceMetric class defined in scikit-learn's library lists all available
        metrics.

    Returns
    -------
    median_min_dist : float
        The median of the distribution of distances between nearest-neighbors.
    """

    data = np.atleast_2d(data)

    nearest_distances = kneighbors_graph(data, 1, mode='distance', metric=metric, include_self=False).data

    median_min_dist = np.median(nearest_distances, overwrite_input=True)

    return round(median_min_dist, 4)


def get_local_densities(data, kernel_mult=2.0, metric='manhattan'):
    """For each sample point of the data-set 'data', estimate a local density in feature
        space by counting the number of neighboring data-points within a particular
        region centered around that sample point.

    Parameters
    ----------
    data : array of shape (n_samples, n_features)
        The data-set, a fraction of whose sample points will be extracted
        by density sampling.

    kernel_mult : float, optional (default = 2.0)
        The kernel multiplier, which determine (in terms of the median of the distribution
        of distances among nearest neighbors) the extent of the regions centered
        around each sample point to consider for the computation of the local density
        associated to that particular sample point.

    metric : string, optional (default = 'manhattan')
        The distance metric used to determine the nearest-neighbor to each data-point.
        The DistanceMetric class defined in scikit-learn's library lists all available
        metrics.

    Returns
    -------
    local_densities : array of shape (n_samples,)
        The i-th entry of this vector corresponds to the local density of the i-th sample
        point in the order of the rows of 'data'.
    """

    data = np.atleast_2d(data)

    assert isinstance(kernel_mult, numbers.Real) and kernel_mult > 0

    kernel_width = kernel_mult * median_min_distance(data, metric)

    N_samples = data.shape[0]

    if 8.0 * get_chunk_size(N_samples, 1) > N_samples:
        A = radius_neighbors_graph(data, kernel_width, mode='connectivity', metric=metric, include_self=True)

        rows, _ = A.nonzero()
        # remove the memmap array usage for now, has an issue...  and why do we need it anyway?
        # with NamedTemporaryFile('w', delete = True, dir = './') as file_name:
        #     fp = np.memmap(file_name, dtype = int, mode = 'w+', shape = rows.shape)
        #     fp[:] = rows[:]
        #     _, counts = np.unique(fp, return_counts = True)
        _, counts = np.unique(rows, return_counts=True)

        local_densities = np.zeros(N_samples, dtype=int)
        for i in range(N_samples):  # CARSTEN
            local_densities[i] = counts[i]
    else:
        local_densities = np.zeros(N_samples, dtype=int)

        chunks_size = get_chunk_size(N_samples, 2)
        for i in range(0, N_samples, chunks_size):
            chunk = data[i:min(i + chunks_size, N_samples)]

            D = pairwise_distances(chunk, data, metric, n_jobs=1)

            D = (D <= kernel_width)

            local_densities[i + np.arange(min(chunks_size, N_samples - i))] = D.sum(axis=1)

    return local_densities


def density_sampling(data, local_densities=None, metric='manhattan',
                     kernel_mult=2.0, outlier_percentile=0.01,
                     target_percentile=0.05, desired_samples=None):
    """The i-th sample point of the data-set 'data' is selected by density sampling
        with a probability given by:

                                      | 0 if outlier_density > LD[i];
        P(keep the i-th data-point) = | 1 if outlier_density <= LD[i] <= target_density;
                                      | target_density / LD[i] if LD[i] > target_density.

        Here 'LD[i]' denotes the local density of the i-th sample point of the data-set,
        whereas 'outlier_density' and 'target_density' are computed as particular percentiles
        of that distribution of local densities.

    Parameters
    ----------
    data : array of shape (n_samples, n_features)
        The data-set, a fraction of whose sample points will be extracted
        by density sampling.

    local_densities : array of shape (n_samples,), optional (default = None)
        The i-th entry of this vector corresponds to the local density of the i-th sample
        point in the order of the rows of 'data'.

    metric : string, optional (default = 'manhattan')
        The distance metric used to determine the nearest-neighbor to each data-point.
        The DistanceMetric class defined in scikit-learn's library lists all available
        metrics.

    kernel_mult : float, optional (default = 2.0)
        The kernel multiplier, which determine (in terms of the median of the distribution
        of distances among nearest neighbors) the extent of the regions centered
        around each sample point to consider for the computation of the local density
        associated to that particular sample point.

    outlier_percentile : float, optional (default = 0.01)
        Specify the outlier density as a percentile of the distribution of local densities.

    target_percentile : float, optional (default = 0.05)
        Specifiy the target density as a percentile of the distribution of local densities.
        Relevant only if 'desired_samples' is left unspecified.

    desired_samples : int, optional (default = None)
        The number of samples to be selected from the whole data-set such that members
        of rare populations and members of more common populations are roughly
        equally represented. To that purpose, a target density is computed that to selects about
        'desired_samples' data-points.

    Returns
    -------
    samples_kept : array of shape (n_selected_samples,)
        If the 'i'-th sample point of 'data' has been selected by a given instance of
        density sampling, number 'i' is featured in the array returned by
        the present function.
    """

    random_state = np.random.RandomState()

    data = np.atleast_2d(data)

    for x in (kernel_mult, outlier_percentile, target_percentile):
        assert isinstance(x, numbers.Real) and x > 0
    for x in (outlier_percentile, target_percentile):
        assert x <= 1.0

    if local_densities is None:
        local_densities = get_local_densities(data, kernel_mult, metric)

    if reduce(operator.mul, local_densities.shape, 1) != max(local_densities.shape):
        raise ValueError("\nERROR: Density_Sampling: density_sampling: problem with "
                         "the dimensions of the vector of local densities provided.\n")
    else:
        local_densities = np.reshape(local_densities, local_densities.size)

    outlier_density = np.percentile(local_densities, outlier_percentile)
    target_density = np.percentile(local_densities, target_percentile)

    samples_kept = np.where(local_densities > outlier_density)[0]
    N_kept = samples_kept.size

    local_densities = local_densities[samples_kept]

    if desired_samples is None:
        probs = np.divide(target_density + 0.0, local_densities)
        ind = np.where(probs > random_state.uniform(size=N_kept))[0]
        samples_kept = samples_kept[ind]
    elif desired_samples <= N_kept:
        sorted_densities = np.sort(local_densities)

        temp = np.reciprocal(sorted_densities[::-1].astype(float))
        cdf = np.cumsum(temp)[::-1]

        target_density = (desired_samples + 0.0) / cdf[0]
        if target_density > sorted_densities[0]:
            temp = desired_samples - np.arange(1.0, N_kept + 1.0)
            possible_targets = np.divide(temp, cdf)

            ind = np.argmax(possible_targets < sorted_densities)
            target_density = possible_targets[ind]

        probs = np.divide(target_density + 0.0, local_densities)
        ind = np.where(probs > random_state.uniform(size=N_kept))[0]
        samples_kept = samples_kept[ind]
    else:
        # print("\nERROR: Density_Sampling: density_sampling: 'desired_samples' has been "
        #       "assigned a value of {desired_samples}, larger than {N_kept}, "
        #       "the number of samples whose local densities are high enough "
        #       "(i.e. excluded are the local densities in the lowest {outlier_percentile} "
        #       "percentile).\n".format(**locals()))
        # exit(1)

        # then just keep them all and don't complain:
        return samples_kept

    return samples_kept


if __name__ == '__main__':
    import doctest
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import datasets
    from sklearn.decomposition import PCA
    from time import sleep


    def plot_PCA(X_reduced, Y, title):
        fig = plt.figure(1, figsize=(10, 8))
        ax = Axes3D(fig, elev=-150, azim=110)

        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=Y,
                   cmap=plt.cm.Paired)

        ax.set_title('First three PCA direction for {title}'.format(**locals()))
        ax.set_xlabel('1st eigenvector')
        ax.w_xaxis.set_ticklabels([])
        ax.set_ylabel('2nd eigenvector')
        ax.w_yaxis.set_ticklabels([])
        ax.set_zlabel('3rd eigenvector')
        ax.w_zaxis.set_ticklabels([])

        plt.show(block=False)
        sleep(3)
        plt.close()


    doctest.testmod()
