import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.offline as py
import cPickle as pickle
save=0
load=1
if save:
    df = pd.read_csv('clusteringMatrix3.csv',index_col=0)
    scaffolds = list(df.axes[0])
    pickle.dump(scaffolds,open('scaffolds.p','wb'))
    data = df.as_matrix()
    pca = PCA(n_components=3)
    pca_transformed = pca.fit_transform(data)
    del data
    np.save('transformed_matrix.npy',pca_transformed)
if load:
    scaffolds = pickle.load(open('scaffolds.p','rb'))
    pca_transformed = np.load('transformed_matrix.npy')

#pd.DataFrame(pca_transformed).to_csv('Transformed_ClusteringMatrix.csv',index=False)
#pca_transformed = pd.read_csv('Transformed_ClusteringMatrix.csv',index_col=False).to_matrix
#
def runKmeans(i):
    global pca_transformed, scaffolds
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



    return inertia

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