import pandas as pd, numpy as np, cPickle as pickle
from sklearn.cluster import FeatureAgglomeration
from sklearn.random_projection import GaussianRandomProjection
from sklearn.decomposition import PCA, KernelPCA, FactorAnalysis, LatentDirichletAllocation, TruncatedSVD,ProjectedGradientNMF, FastICA
from sklearn.manifold import TSNE, SpectralEmbedding
import plotly.graph_objs as go
import plotly.offline as py
save=1
dimensionalityReducers = {'pca':PCA(n_components=3),'kpca':KernelPCA(n_components=3),'feature':FeatureAgglomeration(n_clusters=3),'gauss':GaussianRandomProjection(n_components=3),'factor':FactorAnalysis(n_components=3),'latent':LatentDirichletAllocation(n_topics=3),'SVD':TruncatedSVD(n_components=3),'proj':ProjectedGradientNMF(n_components=3),'ica':FastICA(n_components=3)}
dimensionalityReducers = {'spectral': SpectralEmbedding(n_components=3),'tsne':TSNE(n_components=3)}
if save:
    #df = pd.read_csv('clusteringMatrix3.csv',index_col=0)
    df = pd.read_feather('clusteringMatrix.feather')
    df = df.set_index(['index'])
    scaffolds = list(df.axes[0])
    data = df.as_matrix()
    #np.save('original_matrix.npy',data)
    del df
    for model in dimensionalityReducers:
        try:
            transformed_data = dimensionalityReducers[model].fit_transform(data)
            np.save('transformed3D'+model+'.npy',transformed_data)
            plots = []
            plots.append(go.Scatter3d(x=transformed_data[:, 0], y=transformed_data[:, 1], z=transformed_data[:, 2], name='Data', mode='markers',
                                      marker=dict(color='b', size=2), text=scaffolds))
            fig = go.Figure(data=plots)
            py.plot(fig, filename=model + 'ReduceTest.html')
        except:
            print 'Error with ' + model
    #pca = KernelPCA(n_components=2)#PCA
    #pca_transformed = pca.fit_transform(data)
    #np.save('transformed_matrix2D_kpca.npy' , pca_transformed)
    """
    if tsne_transform:
        transform = 'tsne'
        tsne = TSNE(n_components=n_subgenomes)
        pca_transformed = tsne.fit_transform(data)
        np.save('transformed_matrix%dD_%s.npy' % (n_subgenomes, transform), pca_transformed)"""
    del data