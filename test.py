#coding:utf-8
"""
Example code for creating and visualizing
cluster of county-based cancer risk data

Note that you must download the file
http://www.codeskulptor.org/#alg_clusters_matplotlib.py
to use the matplotlib version of this code
"""

# Flavor of Python - desktop or CodeSkulptor
DESKTOP = True

import math
import random
import urllib2
import alg_cluster

# conditional imports
if DESKTOP:
    import alg_project3_solution      # desktop project solution
    import alg_clusters_matplotlib
else:
    #import userXX_XXXXXXXX as alg_project3_solution   # CodeSkulptor project solution
    import alg_clusters_simplegui
    import codeskulptor
    codeskulptor.set_timeout(30)


###################################################
# Code to load data tables

# URLs for cancer risk data tables of various sizes
# Numbers indicate number of counties in data table

# DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
# DATA_3108_URL = DIRECTORY + "data_clustering/unifiedCancerData_3108.csv"
# DATA_896_URL = DIRECTORY + "data_clustering/unifiedCancerData_896.csv"
# DATA_290_URL = DIRECTORY + "data_clustering/unifiedCancerData_290.csv"
# DATA_111_URL = DIRECTORY + "data_clustering/unifiedCancerData_111.csv"

DATA_3108_URL = "(DATA_3108_URL)unifiedCancerData_3108.csv"
DATA_896_URL = "(DATA_896_URL)unifiedCancerData_896.csv"
DATA_290_URL = "(DATA_290_URL)unifiedCancerData_290.csv"
DATA_111_URL = "(DATA_111_URL)unifiedCancerData_111.csv"

def load_data_table(data_url):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """

    #data_file = urllib2.urlopen(data_url)
    data_file=open('dat/'+data_url)
    data = data_file.read()
    data_file.close()
    data_lines = data.split('\n')
    print "Loaded", len(data_lines), "data points"
    data_tokens = [line.split(',') for line in data_lines]
    return [[tokens[0], float(tokens[1]), float(tokens[2]), int(tokens[3]), float(tokens[4])] 
            for tokens in data_tokens]


############################################################
# Code to create sequential clustering
# Create alphabetical clusters for county data

def sequential_clustering(singleton_list, num_clusters):
    """
    Take a data table and create a list of clusters
    by partitioning the table into clusters based on its ordering
    
    Note that method may return num_clusters or num_clusters + 1 final clusters
    """
    
    cluster_list = []
    cluster_idx = 0
    total_clusters = len(singleton_list)
    cluster_size = float(total_clusters)  / num_clusters
    
    for cluster_idx in range(len(singleton_list)):
        new_cluster = singleton_list[cluster_idx]
        if math.floor(cluster_idx / cluster_size) != \
           math.floor((cluster_idx - 1) / cluster_size):
            cluster_list.append(new_cluster)
        else:
            cluster_list[-1] = cluster_list[-1].merge_clusters(new_cluster)
            
    return cluster_list
                

#####################################################################
# Code to load cancer data, compute a clustering and 
# visualize the results


def run_example():
    """
    Load a data table, compute a list of clusters and 
    plot a list of clusters

    Set DESKTOP = True/False to use either matplotlib or simplegui
    """
    # data_table = load_data_table(DATA_3108_URL)

    print 'in run_example'

    k_n=[]
    h_n=[]
    for x in range(6,21):
        print '------>:',x,'<-----\n'
        # kmeans
        data_table=load_data_table(DATA_111_URL)
        singleton_list=[]
        for line in data_table:
            singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
        cluster_list_k = alg_project3_solution.kmeans_clustering(singleton_list, x, 5)
        kmeans=reduce(lambda x,y:x+y,map(lambda x:x.cluster_error(data_table),cluster_list_k))
        k_n.append(kmeans)

        #hierarchical
        data_table=load_data_table(DATA_111_URL)
        singleton_list=[]
        for line in data_table:
            singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
        cluster_list_h = alg_project3_solution.hierarchical_clustering(singleton_list, x)
        hierarchical=reduce(lambda x,y:x+y,map(lambda x:x.cluster_error(data_table),cluster_list_h))
        h_n.append(hierarchical)
    print 'kmean:',k_n
    print 'hierarchical:',h_n







    # draw the clusters using matplotlib or simplegui
    # if DESKTOP:
    #     # alg_clusters_matplotlib.plot_clusters(data_table, cluster_list, False)
    #     alg_clusters_matplotlib.plot_clusters(data_table, cluster_list_h, True)  #add cluster centers
    # else:
    #     alg_clusters_simplegui.PlotClusters(data_table, cluster_list_h)   # use toggle in GUI to add cluster centers
def get_time():
    import time
    return time.clock()


run_example()





    





  
        






        




