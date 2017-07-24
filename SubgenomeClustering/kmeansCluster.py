import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from multiprocessing import Pool
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import random
import plotly.offline as py
df = pd.read_csv('clusteringMatrix3.csv',index_col=0)
scaffolds = list(df.axes[0])

data = df.as_matrix()
pca = PCA(n_components=3)
pca_transformed = pca.fit_transform(data)
def runKmeans(i):
    global pca_transformed, scaffolds

    kmeans = KMeans(n_clusters=i)
    kmeans.fit(pca_transformed)
    labels = kmeans.labels_
    centers = kmeans.cluster_centers_
    print centers, kmeans.inertia_
    with open('kmeans_nclusters_%d.txt'%i,'w') as f:
        f.write('\n'.join('%s\t%d\tdistance=%f'%(scaffolds[j],labels[j],np.linalg.norm(pca_transformed[j,:]-centers[labels[j],:])) for j in range(len(scaffolds))))
    plots = []
    N = len(labels)
    c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
    for key in labels:
        if list(np.array(scaffolds)[labels == key]):
            for scaffold in np.array(scaffolds)[labels == key]:
                x = pca_transformed[labels == key, 0]
                y = pca_transformed[labels == key, 1]
                z = pca_transformed[labels == key, 2]
                plots.append(go.Scatter3d(x=x,y=y,z=z,name = 'Cluster %d'%key,mode='markers',marker=dict(color=c[key],size=2),text=scaffold))
    x = centers[:, 0]
    y = centers[:, 1]
    z = centers[:, 2]
    plots.append(
        go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(color='purple', symbol='circle', size=12), opacity=0.4,
                     name='Centroids'))
    fig = go.Figure(data=plots)
    py.plot(fig, filename='kmeansOut%d.html' % i)



    return kmeans.inertia_

if __name__ == '__main__':
    p = Pool(9)
    inertias = p.map(runKmeans,range(2,11))
    p.close()
    p.join()


fig = plt.figure()
plt.plot(range(2,11),inertias)
plt.xlabel('Number of Clusters')
plt.ylabel('Inertia of Model')
plt.title('Number of Clusters Versus Inertia')
plt.savefig('SubgenomeClusterModelInertia.png')