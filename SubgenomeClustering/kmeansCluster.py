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
save=0
load=1
n_subgenomes = 3
pca_transform, tsne_transform = 1,0
HDBSCAN,kmeans,gaussian, agglomerate = 0,1,0,1
clusterSizes = range(15,85,10)

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
        pca = KernelPCA(n_components=n_subgenomes)#PCA
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
    pca_transformed = np.load('transformed_matrix%dD_%s.npy'%(n_subgenomes,transform))#.npy')#

#pd.DataFrame(pca_transformed).to_csv('Transformed_ClusteringMatrix.csv',index=False)
#pca_transformed = pd.read_csv('Transformed_ClusteringMatrix.csv',index_col=False).to_matrix
#

# reexpress Data
#z = pca_transformed[:,-1]
#pca_transformed[z>=0,-1] = np.sin(pca_transformed[z>=0,-1])
#pca_transformed[z<0,-1] = np.cos(pca_transformed[z<0,-1])
#pca_transformed = stats.boxcox(pca_transformed- np.amin(pca_transformed,axis=0)+.0001)
# np.concatenate((pca_transformed[:,0],np.vectorize(lambda x: x**(-2))(pca_transformed[:,1]),pca_transformed[:,2]),axis=1)#np.vectorize(lambda x: x**(-2))(pca_transformed)#np.arctan(pca_transformed)# - np.amin(pca_transformed,axis=0)+.0001)

#metrics = ['braycurtis','canberra','chebyshev','cityblock','dice','euclidean','minkowski','p','pyfunc']
#metrics=['braycurtis']
pca_transformed = MinMaxScaler().fit_transform(StandardScaler().fit_transform(pca_transformed))#np.gradient(),axis=1)
print pca_transformed

def runHDBSCAN(clusterSize,metric='euclidean',pca_transformed=pca_transformed,scaffolds=scaffolds):
    #clusterer = hdbscan.HDBSCAN(min_cluster_size=clusterSize,metric = metric,alpha=0.6)#,cluster_selection_method='leaf')#min_samples=5,
    clusterer = DBSCAN(eps=1.2,min_samples=clusterSize)
    clusterer.fit(pca_transformed)
    labels = clusterer.labels_#.fit_predict(pca_transformed)
    #labels = clusterer.labels_
    scores = np.ones(len(scaffolds))#clusterer.probabilities_
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
            if n_subgenomes > 2:
                z = pca_transformed[labels == key, 2]
                plots.append(go.Scatter3d(x=x,y=y,z=z ,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key+1],size=2),text=cluster_scaffolds))#,text=scaffold [idx]
            else:
                plots.append(go.Scatter(x=x,y=y ,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key+1],size=2),text=cluster_scaffolds))#,text=scaffold [idx]

            #Scatter3d with z=z
            fig = go.Figure(data=plots)
            py.plot(fig, filename='hdbscanOut_%s_cluster_%d.html'%(metric,clusterSize))
    return [clusterSize,list(clusterer.labels_).count(-1)/float(len(clusterer.labels_)), len(set(clusterer.labels_))]

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
            #for idx, scaffold in enumerate(cluster_scaffolds):
            #   #bar.update(idx)
            x = pca_transformed[labels == key, 0]
            y = pca_transformed[labels == key, 1]
            if n_subgenomes > 2:
                z = pca_transformed[labels == key, 2]

                plots.append(go.Scatter3d(x=x,y=y,z=z,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key],size=2),text=cluster_scaffolds))#,text=scaffold [idx]
            else:
                plots.append(go.Scatter(x=x,y=y,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key],size=2),text=cluster_scaffolds))#,text=scaffold [idx]

    x = centers[:, 0]
    y = centers[:, 1]
    if n_subgenomes > 2:
        z = centers[:, 2]
        plots.append(
            go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color='purple', symbol='circle', size=12), opacity=0.4,
                         name='Centroids'))
    else:
        plots.append(
            go.Scatter(x=x, y=y, mode='markers', marker=dict(color='purple', symbol='circle', size=12),
                         opacity=0.4,
                         name='Centroids'))
    fig = go.Figure(data=plots)
    py.plot(fig, filename='kmeansOut%d.html' % i)

    return inertia


def findAllLowest(d):
    """search through this dictionary and return all values at the bottom level"""
    for d_list in d:
        print d_list
        findAllLowest(d_list[0])


def runAgglomerate(i,pca_transformed=pca_transformed,scaffolds=scaffolds):
    n = 200

    kmeans = KMeans(n_clusters=n)
    kmeans.fit(pca_transformed)
    inertia = kmeans.inertia_
    labels = kmeans.labels_
    centers = kmeans.cluster_centers_
    origData = defaultdict(list)
    for j in range(len(labels)):
        origData[labels[j]].append(j)
    del kmeans
    print centers, inertia
    """
    origData = defaultdict(list)
    centers = pca_transformed
    for j in range(len(scaffolds)):
        origData[j].append(j)"""
    for m in [i]:#range(n-1,i+2,-10)+[i]:#range(50-1,i+2,-2) + [i]:
        print m, int(m*2/3)
        connectivity = kneighbors_graph(centers, 30, include_self=False)
        agg = AgglomerativeClustering(n_clusters=m,connectivity=connectivity)
        agg.fit(centers)
        labels = agg.labels_
        center_dict = defaultdict(list)
        higherLevelData = defaultdict(list)
        for j in range(len(labels)):
            higherLevelData[labels[j]] += origData[j]
            center_dict[labels[j]].append(centers[labels[j]])
        origData = {k:list(set(higherLevelData[k])) for k in higherLevelData}
        print center_dict
        centers = np.array([np.mean(np.array(center_dict[key]),axis=0) for key in center_dict])
        print centers, 'hi'
    print origData
    labels_final = np.zeros(len(scaffolds))
    for k in origData:
        for val in origData[k]:
            labels_final[val] = k
    labels = np.vectorize(lambda x: int(x))(labels_final)
    print pca_transformed,centers,labels
    inertia = np.sum(np.array([np.linalg.norm(pca_transformed[j, :] - centers[labels[j], :]) for j in range(len(scaffolds))]) ** 2) ** (1 / 2)
    with open('agglomerative_nclusters_%d.txt'%i,'w') as f:
        f.write('inertia=%s\nclusterCenters:\n%s\n'%tuple(map(str,[inertia,centers]))+'\n'.join('%s\t%d\tdistance=%f'%(scaffolds[j],labels[j],np.linalg.norm(pca_transformed[j,:]-centers[labels[j],:])) for j in range(len(scaffolds))))

    """
    connectivity = kneighbors_graph(centers, 2, include_self=False)
    agg = AgglomerativeClustering(n_clusters=i,connectivity=connectivity)
    agg.fit(centers)
    labels = agg.labels_
    final_labels = defaultdict(list)

    higherLevelData = defaultdict(list)
    for j in range(len(labels)):
        higherLevelData[labels[j]].append(origData[j])
    origData = {k: list(set(higherLevelData[k])) for k in higherLevelData}
    print higherLevelData

    for k in higherLevelData:
        final_labels[k].append(findAllLowest(higherLevelData[k]))
        """


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
            #for idx, scaffold in enumerate(cluster_scaffolds):
            #   #bar.update(idx)
            x = pca_transformed[labels == key, 0]
            y = pca_transformed[labels == key, 1]
            if n_subgenomes > 2:
                z = pca_transformed[labels == key, 2]

                plots.append(go.Scatter3d(x=x,y=y,z=z,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key],size=2),text=cluster_scaffolds))#,text=scaffold [idx]
            else:
                plots.append(go.Scatter(x=x,y=y,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key],size=2),text=cluster_scaffolds))#,text=scaffold [idx]

    x = centers[:, 0]
    y = centers[:, 1]
    if n_subgenomes > 2:
        z = centers[:, 2]
        plots.append(
            go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color='purple', symbol='circle', size=12), opacity=0.4,
                         name='Centroids'))
    else:
        plots.append(
            go.Scatter(x=x, y=y, mode='markers', marker=dict(color='purple', symbol='circle', size=12),
                         opacity=0.4,
                         name='Centroids'))
    fig = go.Figure(data=plots)
    py.plot(fig, filename='AgglomerativeOut%d.html' % i)

    return inertia

def runGaussianModel(i,data,pca_transformed=pca_transformed,scaffolds=scaffolds):
    kmeans = mixture.BayesianGaussianMixture(
        n_components=i, covariance_type='full', weight_concentration_prior=1e+2,
        weight_concentration_prior_type='dirichlet_process',
        mean_precision_prior=1e-2, covariance_prior=1e0 * np.eye(2),
        init_params="kmeans", max_iter=100, random_state=2)#kmeans

    kmeans.fit(data)

    labels = kmeans.predict(data)
    centers = kmeans.means_
    inertia = np.sum(np.array([np.linalg.norm(data[j, :] - centers[labels[j], :]) for j in range(len(scaffolds))])**2)**(1/2)
    del kmeans
    print centers, inertia
    with open('gaussian_nclusters_%d.txt' % i, 'w') as f:
        f.write('inertia=%s\nclusterCenters:\n%s\n' % tuple(map(str, [inertia, centers])) + '\n'.join(
            '%s\t%d\tdistance=%f' % (
            scaffolds[j], labels[j], np.linalg.norm(data[j, :] - centers[labels[j], :])) for j in
            range(len(scaffolds))))
    plots = []
    N = len(set(labels))
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N+1)]
    print c
    for key in set(labels):
        print key
        cluster_scaffolds = np.array(scaffolds)[labels == key]
        if list(cluster_scaffolds):
            # bar = progressbar.ProgressBar(max_value=len(cluster_scaffolds),
            #                              widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (',
            #                                       progressbar.ETA(), ') ', ])
            # for idx, scaffold in enumerate(cluster_scaffolds):
            #   #bar.update(idx)
            x = pca_transformed[labels == key, 0]
            y = pca_transformed[labels == key, 1]
            if n_subgenomes > 2:
                z = pca_transformed[labels == key, 2]

                plots.append(go.Scatter3d(x=x, y=y, z=z, name='Cluster %d' % key, mode='markers',
                                          marker=dict(color=c[key], size=2),
                                          text=cluster_scaffolds))  # ,text=scaffold [idx]
            else:
                plots.append(
                    go.Scatter(x=x, y=y, name='Cluster %d' % key, mode='markers', marker=dict(color=c[key], size=2),
                               text=cluster_scaffolds))  # ,text=scaffold [idx]
    """
    x = centers[:, 0]
    y = centers[:, 1]
    if n_subgenomes > 2:
        z = centers[:, 2]
        plots.append(
            go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color='purple', symbol='circle', size=12),
                         opacity=0.4,
                         name='Centroids'))
    else:
        plots.append(
            go.Scatter(x=x, y=y, mode='markers', marker=dict(color='purple', symbol='circle', size=12),
                       opacity=0.4,
                       name='Centroids'))
    """
    fig = go.Figure(data=plots)
    py.plot(fig, filename='gaussianOut%d.html' % i)

    return inertia




if HDBSCAN:
    outData = np.array([runHDBSCAN(clusterSize) for clusterSize in clusterSizes])
    print outData
    fig = plt.figure()
    plt.subplot(211)
    plt.plot(outData[:,0], outData[:,1])
    plt.xlabel('Min Cluster Size')
    plt.ylabel('Proportion Outliers')
    plt.title('Proportion Outliers and Number Clusters versus Min Cluster Size')
    plt.subplot(212)
    plt.xlabel('Min Cluster Size')
    plt.ylabel('Number Clusters')
    plt.plot(outData[:, 0], outData[:, 2])

    plt.savefig('SubgenomeClusterHDBSCAN.png')


if agglomerate:
    #if __name__ == '__main__':
    #    p = Pool()
    inertias = []
    for i in range(2,8):
        inertias.append(runAgglomerate(i))
    #    p.close()
    #    p.join()


    fig = plt.figure()
    plt.plot(range(2,8),inertias)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Inertia of Model')
    plt.title('Number of Clusters Versus Inertia')
    plt.savefig('SubgenomeClusterAgglomerativeModelInertia.png')

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
    plt.savefig('SubgenomeClusterKMeansModelInertia.png')

if gaussian:
    # if __name__ == '__main__':
    #    p = Pool()
    inertias = []
    if n_subgenomes > 2:
        pca = PCA(n_components=2,whiten=True)
        data = pca.fit_transform(pca_transformed)
        data = MinMaxScaler().fit_transform(data)
        #data = data - np.min(data)
        #data = np.array(zip(np.arctan(data[:,1]/data[:,0]),np.zeros(len(data[:,0]))))
    else:
        data = pca_transformed
        data = MinMaxScaler().fit_transform(data)
        data = np.array(zip(np.arctan(data[:, 1] / data[:, 0]), np.zeros(len(data[:, 0]))))

    for i in range(2, 8):
        inertias.append(runGaussianModel(i,data))
        #    p.close()
        #    p.join()


    fig = plt.figure()
    plt.plot(range(2,8),inertias)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Inertia of Model')
    plt.title('Number of Clusters Versus Inertia')
    plt.savefig('SubgenomeClusterGaussianModelInertia.png')