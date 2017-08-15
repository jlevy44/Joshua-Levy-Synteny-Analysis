from sklearn.neural_network import BernoulliRBM as RBM
from sklearn import linear_model, datasets, metrics
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as py
from matplotlib import pyplot as plt
import cPickle as pickle

scaffolds = pickle.load(open('scaffolds.p','rb'))
pca_transformed = np.load('transformed_matrix3D_pca.npy')#[np.vectorize(lambda x: x.startswith('BhS'))(scaffolds)]

logistic = linear_model.LogisticRegression()
rbm = RBM(random_state=0,verbose=True,n_components=3)
rbm.fit(pca_transformed)
print rbm.intercept_hidden_
print rbm.intercept_visible_
print rbm.components_
print pca_transformed
a=rbm.transform(pca_transformed)
plots = []
plots.append(go.Scatter3d(x=a[:,0],y=a[:,1],z=a[:,2] ,name = 'Cluster %d'%1,mode='markers',marker=dict(color='b',size=2),text=scaffolds))
#plt.figure()
#plt.scatter(a[:,0],a[:,1])
#plt.savefig('out.png')
fig = go.Figure(data=plots)
py.plot(fig, filename='deepLearnTest.html')
print np.array(scaffolds) + np.vectorize(lambda x: str(x))(rbm.score_samples(pca_transformed))

"""
classifier = Pipeline(steps=[('rbm',rbm),('logistic',logistic)])
rbm.learning_rate = 0.06
rbm.n_iter = 20

rbm.n_components = 100
logistic.C = 6000.0

classifier.fit(pca_transformed)

logistic_classifier = linear_model.LogisticRegression(C=100.0)
logistic_classifier.fit(pca_transformed)

print("Logistic regression using RBM features:\n%s\n" % (
    metrics.classification_report(
        classifier.predict(pca_transformed))))

print("Logistic regression using raw pixel features:\n%s\n" % (
    metrics.classification_report(
        logistic_classifier.predict(pca_transformed))))"""