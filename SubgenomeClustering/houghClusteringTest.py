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
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from scipy import stats
import cv2
from gmr import GMM, plot_error_ellipses

def hough_line(centers):
  # Rho and Theta ranges
  thetas = np.deg2rad(np.arange(-90.0, 90.0))
  width, height = (np.max(centers[:,0])-np.min(centers[:,0]),np.max(centers[:,1])-np.min(centers[:,1]))
  diag_len = int(np.ceil(np.sqrt(width * width + height * height)))   # max_dist
  rhos = np.linspace(-diag_len, diag_len, diag_len * 2.0)

  # Cache some resuable values
  cos_t = np.cos(thetas)
  sin_t = np.sin(thetas)
  num_thetas = len(thetas)

  # Hough accumulator array of theta vs rho
  #print diag_len,num_thetas
  accumulator = np.zeros((2 * diag_len, num_thetas), dtype=np.uint64)
  #y_idxs, x_idxs = np.nonzero(img)  # (row, col) indexes to edges

  # Vote in the hough accumulator
  for i in range(len(centers)):
    #x = x_idxs[i]
    #y = y_idxs[i]
    x = centers[i,0]
    y = centers[i, 1]

    for t_idx in range(num_thetas):
      # Calculate rho. diag_len is added for a positive index
      rho = int(round(x * cos_t[t_idx] + y * sin_t[t_idx]) + diag_len)

      accumulator[rho, t_idx] += 1


  return accumulator, thetas, rhos
# code adopted  from online help: alyssaq

def show_hough_line(pca_transformed, thetas, rhos, accumulator):
  fig, ax = plt.subplots(1, 2, figsize=(20, 10))

  #ax[0].imshow(img, cmap=plt.cm.gray)
  ax[0].scatter(pca_transformed[:,0],pca_transformed[:,1],s=5,color='b')
  ax[0].set_title('Input image')
  ax[0].axis([min(pca_transformed[:,0]),max(pca_transformed[:,0]),min(pca_transformed[:,1]),max(pca_transformed[:,1])])

  ax[1].imshow(
    accumulator, cmap='jet',
    extent=[np.rad2deg(thetas[-1]), np.rad2deg(thetas[0]), rhos[-1], rhos[0]])
  ax[1].set_aspect('equal', adjustable='box')
  ax[1].set_title('Hough transform')
  ax[1].set_xlabel('Angles (degrees)')
  ax[1].set_ylabel('Distance (pixels)')
  ax[1].axis('image')

  #plt.axis('off')
  plt.savefig('output.png', bbox_inches='tight')
  plt.show()

scaffolds = pickle.load(open('scaffolds.p','rb'))
pca_transformed = np.load('transformed_matrix2D_pca.npy')#[np.vectorize(lambda x: x.startswith('BhS'))(scaffolds)]#MinMaxScaler().fit_transform(StandardScaler().fit_transform(
scaffolds = pickle.load(open('scaffolds.p','rb'))
connectivity = kneighbors_graph(pca_transformed, 30, include_self=False)
kmeans = AgglomerativeClustering(n_clusters=50,connectivity=connectivity)#KMeans(n_clusters=30)#AgglomerativeClustering(n_clusters=50,connectivity=connectivity)
#kmeans = hdbscan.HDBSCAN()
kmeans.fit(pca_transformed)

#inertia = kmeans.inertia_
center_dict = defaultdict(list)
labels = kmeans.labels_
for j in range(len(labels)):
    center_dict[labels[j]].append(pca_transformed[labels[j]])
centers = np.array([np.mean(np.array(center_dict[key]),axis=0) for key in center_dict])#kmeans.cluster_centers_

edges = centers
"""
lines = cv2.HoughLines(edges,1,np.pi/180,200)
"""
accumulator, thetas, rhos = hough_line(centers)

neighborhood_size = 5
threshold = 70
data_max = filters.maximum_filter(accumulator, neighborhood_size)
maxima = (accumulator == data_max)
data_min = filters.minimum_filter(accumulator, neighborhood_size)
diff = ((data_max - data_min) > threshold)
maxima[diff == 0] = 0
labeled, num_objects = ndimage.label(maxima)
slices = ndimage.find_objects(labeled)
x, y = [], []
for dy,dx in slices:
    x_center = (dx.start + dx.stop - 1)/2
    x.append(x_center)
    y_center = (dy.start + dy.stop - 1)/2
    y.append(y_center)
#print accumulator, rhos, thetas
show_hough_line(centers,thetas,rhos,accumulator)
#lines = cv2.HoughLines(centers,1,np.pi/180,200)
fig = plt.figure()

plt.scatter(pca_transformed[:,0],pca_transformed[:,1],s=1,color='b')
plt.scatter(centers[:,0],centers[:,1],s=50,color='g')
#print np.concatenate((rhos.T,thetas.T),axis=1)
print x,y
for i in range(len(x)):#rho,theta in np.concatenate((rhos,thetas)):
    theta = thetas[x[i]]
    rho = rhos[y[i]]
    a = np.cos(theta)
    b = np.sin(theta)
    x0 = a*rho
    y0 = b*rho
    x1 = int(x0 + 1000*(-b))
    y1 = int(y0 + 1000*(a))
    x2 = int(x0 - 1000*(-b))
    y2 = int(y0 - 1000*(a))
    plt.plot([x0,x1],[y0,y1],color='r')
"""for rho,theta in lines[0]:
    a = np.cos(theta)
    b = np.sin(theta)
    x0 = a*rho
    y0 = b*rho
    x1 = int(x0 + 1000*(-b))
    y1 = int(y0 + 1000*(a))
    x2 = int(x0 - 1000*(-b))
    y2 = int(y0 - 1000*(a))
    plt.plot([x0, x1], [y0, y1], color='r')
    #cv2.line(img,(x1,y1),(x2,y2),(0,0,255),2)"""
plt.axis([-1,2,-1,10])

plt.savefig('outTest.png')

gmm = GMM(n_components=2)
kmeans = KMeans(n_clusters=50)
kmeans.fit(MinMaxScaler().fit_transform(pca_transformed))
centers = kmeans.cluster_centers_
gmm.from_samples(centers)#MinMaxScaler().fit_transform(pca_transformed)


cond = gmm.condition(np.array([0]), np.array([1.0]))

plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.title("Gaussian Mixture Model")
plt.xlim((-10, 10))
plt.ylim((-10, 10))
plot_error_ellipses(plt.gca(), gmm, colors=["r", "g", "b"])
plt.scatter(pca_transformed[:, 0], pca_transformed[:, 1])

plt.subplot(1, 3, 2)
plt.title("Probability Density and Samples")
plt.xlim((-10, 10))
plt.ylim((-10, 10))
x, y = np.meshgrid(np.linspace(-10, 10, 100), np.linspace(-10, 10, 100))
X_test = np.vstack((x.ravel(), y.ravel())).T
p = gmm.to_probability_density(X_test)
p = p.reshape(*x.shape)
plt.contourf(x, y, p)
X_sampled = gmm.sample(100)
plt.scatter(X_sampled[:, 0], X_sampled[:, 1], c="r")

plt.subplot(1, 3, 3)
plt.title("Conditional PDF $p(y | x = 1)$")
X_test = np.linspace(-10, 10, 100)
plt.plot(X_test, cond.to_probability_density(X_test[:, np.newaxis]))

plt.savefig('GMM.png')


"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
pca_transformed = np.load('transformed_matrix2D_pca.npy')
fig = plt.figure()
plt.scatter(pca_transformed[:,0],pca_transformed[:,1],s=1,color='b')
plt.show()
plt.savefig('out2.png')
exit()



"""