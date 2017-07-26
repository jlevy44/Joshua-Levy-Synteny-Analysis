import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.offline as py
import cPickle as pickle
import hdbscan
from scipy import stats
#import progressbar
save=0
load=1
n_subgenomes = 3
pca_transform, tsne_transform = 1,0
HDBSCAN,kmeans = 1,1
if pca_transform:
    transform = 'pca'
else:
    transform= 'tsne'
if save:
    df = pd.read_csv('clusteringMatrix3.csv',index_col=0)
    scaffolds = list(df.axes[0])
    pickle.dump(scaffolds,open('scaffolds.p','wb'))
    data = df.as_matrix()
    del df
    if pca_transform:
        transform = 'pca'
        pca = PCA(n_components=n_subgenomes)
        pca_transformed = pca.fit_transform(data)
        np.save('transformed_matrix%dD_%s.npy' %(n_subgenomes,transform), pca_transformed)
    if tsne_transform:
        transform = 'tsne'
        tsne = TSNE(n_components=n_subgenomes)
        pca_transformed = tsne.fit_transform(data)
        np.save('transformed_matrix%dD_%s.npy' % (n_subgenomes, transform), pca_transformed)
    del data

if load:
    scaffolds = pickle.load(open('scaffolds.p','rb'))
    pca_transformed = np.load('transformed_matrix.npy')#%dD_%s.npy'%(n_subgenomes,transform))

#pd.DataFrame(pca_transformed).to_csv('Transformed_ClusteringMatrix.csv',index=False)
#pca_transformed = pd.read_csv('Transformed_ClusteringMatrix.csv',index_col=False).to_matrix
#

# reexpress Data
#z = pca_transformed[:,-1]
#pca_transformed[z>=0,-1] = np.sin(pca_transformed[z>=0,-1])
#pca_transformed[z<0,-1] = np.cos(pca_transformed[z<0,-1])
#pca_transformed = stats.boxcox(pca_transformed- np.amin(pca_transformed,axis=0)+.0001)
# np.concatenate((pca_transformed[:,0],np.vectorize(lambda x: x**(-2))(pca_transformed[:,1]),pca_transformed[:,2]),axis=1)#np.vectorize(lambda x: x**(-2))(pca_transformed)#np.arctan(pca_transformed)# - np.amin(pca_transformed,axis=0)+.0001)

metrics = ['braycurtis','canberra','chebyshev','cityblock','dice','euclidean','minkowski','p','pyfunc']
metrics=['braycurtis']


def runHDBSCAN(metric,pca_transformed=pca_transformed,scaffolds=scaffolds):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=10,min_samples=0,metric = metric)#,cluster_selection_method='leaf')
    clusterer.fit(pca_transformed)
    labels = clusterer.labels_
    scores = clusterer.probabilities_
    plots = []
    N = len(set(labels))
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
    with open('hdbscan_output.txt','w') as f:
        f.write('\n'.join('%s\t%d\tscore=%f'%(scaffolds[j],labels[j],scores[j]) for j in np.arange(len(scaffolds))))
    text = np.vectorize(lambda i: "%s score=%s"%(scaffolds[i],str(scores[i])))(np.arange(len(scaffolds)))
    for key in set(labels):
        cluster_scaffolds = np.array(text)[labels == key]
        if list(cluster_scaffolds):
            x = pca_transformed[labels == key, 0]
            y = pca_transformed[labels == key, 1]
            z = pca_transformed[labels == key, 2]
            plots.append(go.Scatter3d(x=x,y=y,z=z ,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key+1],size=2),text=cluster_scaffolds))#,text=scaffold [idx]
            #Scatter3d with z=z
            fig = go.Figure(data=plots)
            py.plot(fig, filename='hdbscanOut_%s.html'%metric)

def runKmeans(i,pca_transformed=pca_transformed,scaffolds=scaffolds):
    kmeans = KMeans(n_clusters=i)
    kmeans.fit(pca_transformed)
    inertia = kmeans.inertia_
    labels = kmeans.labels_
    centers = kmeans.cluster_centers_
    del kmeans
    print centers, inertia
    with open('kmeans_nclusters_%d.txt'%i,'w') as f:
        f.write('inertia=%s\nclusterCenters:\n%s\n'%tuple(map(str,[inertia,centers]))+'\n'.join('%s\t%d\tdistance=%f'%(scaffolds[j],labels[j],np.linalg.norm(pca_transformed[j,:]-centers[labels[j],:])) for j in range(len(scaffolds))))
    plots = []
    N = len(set(labels))
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
    for key in set(labels):
        print key
        cluster_scaffolds = np.array(scaffolds)[labels == key]
        if list(cluster_scaffolds):
            #bar = progressbar.ProgressBar(max_value=len(cluster_scaffolds),
            #                              widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (',
            #                                       progressbar.ETA(), ') ', ])
            x = pca_transformed[labels == key, 0]
            y = pca_transformed[labels == key, 1]
            z = pca_transformed[labels == key, 2]
            #for idx, scaffold in enumerate(cluster_scaffolds):
            #   #bar.update(idx)
            plots.append(go.Scatter3d(x=x,y=y,z=z,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key],size=2),text=cluster_scaffolds))#,text=scaffold [idx]
    x = centers[:, 0]
    y = centers[:, 1]
    z = centers[:, 2]
    plots.append(
        go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color='purple', symbol='circle', size=12), opacity=0.4,
                     name='Centroids'))
    fig = go.Figure(data=plots)
    py.plot(fig, filename='kmeansOut%d.html' % i)

    return inertia



if HDBSCAN:
    for metric in metrics:
        runHDBSCAN(metric)

if kmeans:
    #if __name__ == '__main__':
    #    p = Pool()
    inertias = []
    for i in range(2,8):
        inertias.append(runKmeans(i))
    #    p.close()
    #    p.join()


    fig = plt.figure()
    plt.plot(range(2,8),inertias)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Inertia of Model')
    plt.title('Number of Clusters Versus Inertia')
    plt.savefig('SubgenomeClusterModelInertia.png')