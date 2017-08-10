from sklearn.neural_network import BernoulliRBM as RBM
from sklearn import linear_model, datasets, metrics
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn import cluster
import numpy as np
import time
import pandas as pd
import cPickle as pickle
import plotly.graph_objs as go
import plotly.offline as py
import os
from multiprocessing import Pool
from collections import defaultdict
np.random.seed(0)


def classify(classifier,testClusterData, testCluster_scaffolds, trainData, trainData_scaffolds, trainLabels, Tname, n_clusters, name):
    class_name, clf = classifier
    try:
        clf.fit(trainData, trainLabels)
        labels = np.array(clf.predict(testClusterData))
        N = len(set(labels))
        c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
        c2 = {list(set(labels))[i]: c[i] for i in range(N)}
        # below is a test
        plots = []
        for key in set(labels):
            cluster_scaffolds = testCluster_scaffolds[labels == key]
            # print set(labels)
            # print cluster_scaffolds
            # print transformed_data
            print class_name
            print labels
            print key
            # print labels == key
            # print transformed_data[labels==key,0]
            # print len(cluster_scaffolds), len(transformed_data[labels==key,0])
            # print c2
            # print key, labels[labels==key], c2[key]
            # print map(len,[testClusterData[labels==key,0],testClusterData[labels==key,1],testClusterData[labels==key,2],dict(color=c2[key],size=2),cluster_scaffolds])
            if list(cluster_scaffolds):
                plots.append(go.Scatter3d(x=testClusterData[labels == key, 0], y=testClusterData[labels == key, 1],
                                          z=testClusterData[labels == key, 2],
                                          name='Cluster %d, %d points' % (key, len(cluster_scaffolds)), mode='markers',
                                          marker=dict(color=c2[key], size=2), text=cluster_scaffolds))
        fig = go.Figure(data=plots)
        py.plot(fig, filename=name + Tname + 'n%d' % n_clusters + class_name + 'ClassifierClusterTest.html')
        return name + Tname + 'n%d' % n_clusters + class_name + 'ClassifierClusterTest.html'
    except:
        return 'Unable to classify using options %s, %s, %d and %s' % (Tname, name, n_clusters, class_name)

def clusterfy(file):
    global scaffolds
    Tname = file.split('3D')[1].split('.npy')[0]
    transformed_data = np.load(file)

    transformed_data = StandardScaler().fit_transform(transformed_data)

    connectivity = kneighbors_graph(transformed_data, n_neighbors=10, include_self=False)
    connectivity = 0.5 * (connectivity + connectivity.T)
    bandwidth = cluster.estimate_bandwidth(transformed_data, quantile=0.3)


    clustering_names = [
        'MiniBatchKMeans', 'AffinityPropagation', 'MeanShift',
        'SpectralClustering', 'Ward', 'AgglomerativeClustering',
        'DBSCAN', 'Birch']
    n_clusters = 3
    ms = cluster.MeanShift(bandwidth=bandwidth, bin_seeding=True)
    two_means = cluster.MiniBatchKMeans(n_clusters=n_clusters)
    ward = cluster.AgglomerativeClustering(n_clusters=n_clusters, linkage='ward',
                                           connectivity=connectivity)
    spectral = cluster.SpectralClustering(n_clusters=n_clusters,
                                          eigen_solver='arpack',
                                          affinity="nearest_neighbors")#
    dbscan = cluster.DBSCAN(eps=.2)
    affinity_propagation = cluster.AffinityPropagation(damping=.9,
                                                       preference=-200)

    average_linkage = cluster.AgglomerativeClustering(
        linkage="average", affinity="cityblock", n_clusters=n_clusters,
        connectivity=connectivity)

    birch = cluster.Birch(n_clusters=n_clusters)


    clustering_algorithms = [
            two_means, affinity_propagation, ms, spectral, ward, average_linkage,
            dbscan, birch]

    for name, algorithm in zip(clustering_names, clustering_algorithms):
        try:
            # predict cluster memberships
            t0 = time.time()
            algorithm.fit(transformed_data)
            t1 = time.time()
            if hasattr(algorithm, 'labels_'):
                y_pred = algorithm.labels_.astype(np.int)
            else:
                y_pred = algorithm.predict(transformed_data)
            N = len(set(y_pred))
            c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, N)]
            # plot
            plots = []
            clusterSize = defaultdict(list)
            for key in set(y_pred):
                #print key
                cluster_scaffolds = np.array(scaffolds)[y_pred == key]
                if list(cluster_scaffolds):
                    clusterSize[key]=len(cluster_scaffolds)
                    if clusterSize[key] == max(clusterSize.values()):
                        testCluster = key
                    plots.append(go.Scatter3d(x=transformed_data[y_pred==key,0],y=transformed_data[y_pred==key,1],z=transformed_data[y_pred==key,2] ,name = 'Cluster %d, %d points'%(key,len(cluster_scaffolds)),mode='markers',marker=dict(color=c[key],size=2),text=cluster_scaffolds))


            if hasattr(algorithm, 'cluster_centers_'):
                centers = algorithm.cluster_centers_
                plots.append(
                    go.Scatter3d(x=centers[:,0], y=centers[:,1], z=centers[:,2], mode='markers', marker=dict(color='purple', symbol='circle', size=12),
                                 opacity=0.4,
                                 name='Centroids'))

            fig = go.Figure(data=plots)
            py.plot(fig, filename=name + Tname + 'n%d'%n_clusters +'ClusterTest.html')

            #print y_pred == testCluster
            testClusterData = transformed_data[y_pred == testCluster]
            testCluster_scaffolds = np.array(scaffolds)[y_pred == testCluster]
            trainData = transformed_data[y_pred != testCluster]
            trainData_scaffolds = np.array(scaffolds)[y_pred != testCluster]
            trainLabels = y_pred[y_pred != testCluster]

            class_names = ["NearestNeighbors", "LinearSVM", "RBFSVM", "GaussianProcess",
                     "DecisionTree", "RandomForest", "NeuralNet", "AdaBoost",
                     "NaiveBayes", "QDA"]

            classifiers = [
                KNeighborsClassifier(3),
                SVC(kernel="linear", C=0.025,verbose=1),
                SVC(gamma=2, C=1,verbose=1),
                GaussianProcessClassifier(1.0 * RBF(1.0), warm_start=True),
                DecisionTreeClassifier(max_depth=5),
                RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1,verbose=1),
                MLPClassifier(alpha=1,verbose=1),
                AdaBoostClassifier(),
                GaussianNB(),
                QuadraticDiscriminantAnalysis()]

            output = [classify(classifier,testClusterData, testCluster_scaffolds, trainData, trainData_scaffolds, trainLabels, Tname, n_clusters, name) for classifier in zip(class_names,classifiers)]
            with open('classifyOut_%s_%s.txt' % (Tname, name), 'w') as f:
                f.write('\n'.join(output))
        except:
            print 'Unable to cluster completely using ' + name + ' for ' + Tname
    return file



files = [files for files in os.listdir('.') if files.startswith('transformed3D') and files.endswith('.npy')]
scaffolds = pickle.load(open('scaffolds.p', 'rb'))
results = []
def log_result(result):
    results.append(result)

if __name__ == '__main__':
    p = Pool(len(files))
    for file in files:
        p.apply_async(clusterfy, args=(file,), callback=log_result)
    p.close()
    p.join()
    print results
    print 'pool closed '

