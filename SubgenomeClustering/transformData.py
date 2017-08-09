import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.decomposition import KernelPCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, SpectralClustering, SpectralBiclustering, DBSCAN, AgglomerativeClustering
from sklearn import mixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.neighbors import kneighbors_graph
import plotly.graph_objs as go
import plotly.offline as py
from sklearn.preprocessing import MinMaxScaler
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
import cPickle as pickle
import hdbscan
from collections import defaultdict
from scipy import stats
#import progressbar
save=1
load=0
n_subgenomes = 3


if save:
    #df = pd.read_csv('clusteringMatrix3.csv',index_col=0)
    df = pd.read_feather('clusteringMatrix.feather')
    df = df.set_index(['index'])
    data = df.as_matrix()
    #np.save('original_matrix.npy',data)
    del df
    pca = KernelPCA(n_components=2)#PCA
    pca_transformed = pca.fit_transform(data)
    np.save('transformed_matrix2D_kpca.npy' , pca_transformed)
    """
    if tsne_transform:
        transform = 'tsne'
        tsne = TSNE(n_components=n_subgenomes)
        pca_transformed = tsne.fit_transform(data)
        np.save('transformed_matrix%dD_%s.npy' % (n_subgenomes, transform), pca_transformed)"""
    del data

