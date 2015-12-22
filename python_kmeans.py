# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 00:11:35 2015

@author: Harry Yang
"""

import sklearn as learn
from sklearn import cluster
import numpy as np
import scipy
import pandas
from python_kmeans_test_haplotype_generator import haplotype_simulator

def similarity(matrix):
    """
    compare similarity of two matrix. 
    input: 2xlen matrix where len = length of haplo
    output: sim. score (percentage)
    """
    sim_score=0
    for i in range(len(matrix[0])):
        if matrix[0][i] != 0 and matrix[0][i]==matrix[1][i]:
            sim_score=sim_score+1
    return sim_score*1.0/len(matrix[0])

def test(repeat):
    sim_avg=0.0;
    num_successful_clustering=0;
    for i in range(repeat):    
        test_matrix=haplotype_simulator(1000,1000,100,2)
        test_center=[]
        test_center.append(test_matrix[0])
        test_center.append(test_matrix[500])
        sim_avg = sim_avg + similarity(test_center)
#        test_bandwidth=learn.cluster.estimate_bandwidth(test_matrix, quantile=0.5)
        ms_vector = learn.cluster.MeanShift(seeds=test_center).fit_predict(test_matrix)
        if ms_vector[0] == ms_vector[251] and ms_vector[523]==ms_vector[809]:
            num_successful_clustering=num_successful_clustering+1
    print "average simliarty = ", sim_avg/repeat, "and # successful clustering = ", num_successful_clustering
if __name__ == "__main__":
#    print "hello world"
#    t= np.zeros((3,3), dtype = np.double)
#    t[0:1,0]=1
#    t[0:1,1]=0
#    t[0:2.2]=1
#    test=test_matrix_40_reads_100_postxt
#    for i in range(1,101):
#        test[0,i]=test[1,i]
#    k=cluster.k_means(test,n_clusters=2)
#    print k[1]
#    
#    for i in range(len(k[1])):
#        print k[1][i]
    test(10000)
#    print cluster.k_means(test_matrix,n_init=1000,n_clusters=2,max_iter=1000,precompute_distances=True)[1]
#    print cluster.AgglomerativeClustering(n_clusters=2, affinity='euclidean').fit_predict(test_matrix)
#    print cluster.spectral_clustering(test_matrix, n_clusters=2)