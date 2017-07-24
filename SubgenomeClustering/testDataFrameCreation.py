from collections import Counter, defaultdict
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import random
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
#import operator
#from mpl_toolkits_Josh.mplot3d import Axes3D
#import imp
#plot3D = imp.load_source('plot3D','/global/u2/j/jlevy/python_modules/mpl_toolkits/mplot3d')
#Axes3D = plot3D.axes3d()
C = lambda a: Counter(a)
#d = defaultdict(list)
#d['a'] = C(['a','r','t','g','d','f','f','f','d','d','s'])
#d['d'] = C(['a','a','t','t','y','r','w','e'])
#d['e'] = C(['e','r','w','e','i','y','q','w','e','d'])
names = ['John','Andy','Joe','Johnson','Smith','Williams']
d = {i:random.choice(names) for i in range(60)}
d2 = {i:C(list(d[i])) for i in range(60)}
df = pd.DataFrame(d2).T.fillna(0.)
"""
     a    d    e    f    g    i    q    r    s    t    w    y
a  1.0  3.0  0.0  3.0  1.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0
d  2.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0  2.0  1.0  1.0
e  0.0  1.0  3.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  2.0  1.0
"""
rows = list(df.axes[0])
columns = list(df.axes[1])
#df.rows.names = d.values()
df.to_csv('test.csv',index=True)
df2 = pd.read_csv('test.csv',index_col=0)
"""
     a    d    e    f    g    i    q    r    s    t    w    y
a  1.0  3.0  0.0  3.0  1.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0
d  2.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0  2.0  1.0  1.0
e  0.0  1.0  3.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  2.0  1.0
"""
df2Mat = df.as_matrix()
"""
array([[ 1.,  3.,  0.,  3.,  1.,  0.,  0.,  1.,  1.,  1.,  0.,  0.],
       [ 2.,  0.,  1.,  0.,  0.,  0.,  0.,  1.,  0.,  2.,  1.,  1.],
       [ 0.,  1.,  3.,  0.,  0.,  1.,  1.,  1.,  0.,  0.,  2.,  1.]])
"""
data = df2Mat
pca = PCA(n_components=3)
#pca.fit(data)
pca_transformed = pca.fit_transform(data)
inertias = []
for j in range(4,11):
    kmeans = KMeans(n_clusters=j)
    kmeans.fit(pca_transformed)
    inertias.append(kmeans.inertia_)
    centers = kmeans.cluster_centers_
    labels = kmeans.labels_
    #fig = plt.figure()
    colors = {0:'red',1:'blue',2:'yellow',3:'black',4:'cyan',5:'green'}
    #ax = fig.add_subplot(111, projection='3d')
    plots = []
    for key in list(set(labels)):
        if list(np.array(d.values())[labels == key]):
            x = pca_transformed[labels == key,0]
            y = pca_transformed[labels == key,1]
            z = pca_transformed[labels == key,2]
            if key <= 5:
                c = colors[key]
            else:
                c = random.randint(0,255)
            plots.append(go.Scatter3d(x=x,y=y,z=z,name =Counter(np.array(d.values())[labels == key]).most_common()[0][0],mode='markers',marker=dict(color=c,size=5)))
        #ax.scatter(x,y,z,c=colors[key])
    x = kmeans.cluster_centers_[:,0]
    y = kmeans.cluster_centers_[:,1]
    z = kmeans.cluster_centers_[:,2]
    plots.append(go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(color='purple',symbol='circle',size=12),opacity=0.6,name='Centroids'))
    fig = go.Figure(data=plots)
    namesT = d.keys()
    for i in range(len(namesT)):
        print namesT[i],labels[i], np.linalg.norm(pca_transformed[i,:]-centers[labels[i],:])
    py.plot(fig,filename='kmeansTestOut%d.html'%j)
fig = plt.figure()
plt.plot(range(2,11),inertias)
plt.xlabel('Number of Clusters')
plt.ylabel('Inertia of Model')
plt.title('Number of Clusters Versus Inertia')
plt.savefig('ClusterModelInertia.png')

#for idx,label in kmeans.labels_:
