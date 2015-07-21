#coding:utf-8
"""
Cluster class for Module 3
"""

import math


class Cluster:
    """
    Class for creating and merging clusters of counties
    """
    
    def __init__(self, fips_codes, horiz_pos, vert_pos, population, risk):
        """
        Create a cluster based the models a set of counties' data
        """
        self._fips_codes = fips_codes
        self._horiz_center = horiz_pos
        self._vert_center = vert_pos
        self._total_population = population
        self._averaged_risk = risk
        
        
    def __repr__(self):
        """
        String representation assuming the module is "alg_cluster".
        """
        rep = "alg_cluster.Cluster("
        rep += str(self._fips_codes) + ", "
        rep += str(self._horiz_center) + ", "
        rep += str(self._vert_center) + ", "
        rep += str(self._total_population) + ", "
        rep += str(self._averaged_risk) + ")"
        return rep


    def fips_codes(self):
        """
        Get the cluster's set of FIPS codes
        """
        return self._fips_codes
    
    def horiz_center(self):
        """
        Get the averged horizontal center of cluster
        """
        return self._horiz_center
    
    def vert_center(self):
        """
        Get the averaged vertical center of the cluster
        """
        return self._vert_center
    
    def total_population(self):
        """
        Get the total population for the cluster
        """
        return self._total_population
    
    def averaged_risk(self):
        """
        Get the averaged risk for the cluster
        """
        return self._averaged_risk
   
        
    def copy(self):
        """
        Return a copy of a cluster
        """
        copy_cluster = Cluster(set(self._fips_codes), self._horiz_center, self._vert_center,
                               self._total_population, self._averaged_risk)
        return copy_cluster


    def distance(self, other_cluster):
        """
        Compute the Euclidean distance between two clusters
        """
        vert_dist = self._vert_center - other_cluster.vert_center()
        horiz_dist = self._horiz_center - other_cluster.horiz_center()
        return math.sqrt(vert_dist ** 2 + horiz_dist ** 2)
        
    def merge_clusters(self, other_cluster):
        """
        Merge one cluster into another
        The merge uses the relatively populations of each
        cluster in computing a new center and risk
        
        Note that this method mutates self
        """
        if len(other_cluster.fips_codes()) == 0:
            return self
        else:
            self._fips_codes.update(set(other_cluster.fips_codes()))
 
            # compute weights for averaging
            self_weight = float(self._total_population)                        
            other_weight = float(other_cluster.total_population())
            self._total_population = self._total_population + other_cluster.total_population()
            self_weight /= self._total_population
            other_weight /= self._total_population
                    
            # update center and risk using weights
            self._vert_center = self_weight * self._vert_center + other_weight * other_cluster.vert_center()
            self._horiz_center = self_weight * self._horiz_center + other_weight * other_cluster.horiz_center()
            self._averaged_risk = self_weight * self._averaged_risk + other_weight * other_cluster.averaged_risk()
            return self

    def cluster_error(self, data_table):
        """
        Input: data_table is the original table of cancer data used in creating the cluster.
        
        Output: The error as the sum of the square of the distance from each county
        in the cluster to the cluster center (weighted by its population)
        """
        # Build hash table to accelerate error computation
        fips_to_line = {}
        for line_idx in range(len(data_table)):
            line = data_table[line_idx]
            fips_to_line[line[0]] = line_idx
        
        # compute error as weighted squared distance from counties to cluster center
        total_error = 0
        counties = self.fips_codes()
        for county in counties:
            line = data_table[fips_to_line[county]]
            singleton_cluster = Cluster(set([line[0]]), line[1], line[2], line[3], line[4])
            singleton_distance = self.distance(singleton_cluster)
            total_error += (singleton_distance ** 2) * singleton_cluster.total_population()
        return total_error
            
        
            

        
    
    
            


def slow_closest_pair(cluster_list):
    'find closest pair in cluster_list'
    # print 'slow_closest_pair::::::::::::::>',cluster_list,'<:::::::::::::::\n'
    closest_pair=[float('inf'),-1,-1]
    for u_item in xrange(len(cluster_list)):

        for v_item in xrange(u_item+1,len(cluster_list)):
            distance=cluster_list[u_item].distance(cluster_list[v_item])
            if closest_pair[0]>distance:
                closest_pair[0]=distance
                closest_pair[1]=u_item
                closest_pair[2]=v_item
    # print 'return ',closest_pair
    return tuple(closest_pair)

def fast_closest_pair(cluster_list):
    'find closest pair in cluster_list ,quick than slow_closest_pair()'

    # print '\n fast_closest_pair:========>',cluster_list,'<=============\n'
    closest_pair=[float('inf'),-1,-1]
    total_number=len(cluster_list)
    templist=list(cluster_list)
    if total_number<=3:
        closest_pair=slow_closest_pair(cluster_list)
    else:
        mid1=total_number/2
        templist.sort(key = lambda templist: templist.horiz_center())

        cluster_left=templist[:mid1]
        cluster_right=templist[mid1:]
        left_pair=fast_closest_pair(cluster_left)
        right_pair=fast_closest_pair(cluster_right)
        #min(left_pair[0],right_pair[0])
        if left_pair[0]<right_pair[0]:
            closest_pair[0]=left_pair[0]
            closest_pair[1]=cluster_list.index(templist[left_pair[1]])
            closest_pair[2]=cluster_list.index(templist[left_pair[2]])
        else:
            closest_pair[0]=right_pair[0]
            # print mid1,cluster_list,right_pair[1]+mid1
            closest_pair[1]=cluster_list.index(templist[right_pair[1]+mid1])
            closest_pair[2]=cluster_list.index(templist[right_pair[2]+mid1])
        # print float((templist[mid1-1].horiz_center()+templist[mid1].horiz_center()))
        mid2=(1.0/2)*float((templist[mid1-1].horiz_center()+templist[mid1].horiz_center()))

        other_point=closest_pair_strip(cluster_list,mid2,closest_pair[0])
        # print 'other_point:????',other_point,'??????????',other_point[0],other_point[1],other_point[2]
        if closest_pair[0]>other_point[0]:
            closest_pair[0]=other_point[0]

            closest_pair[1]=other_point[1]
            closest_pair[2]=other_point[2]
    # print 'return:',closest_pair
    return tuple(closest_pair)
def closest_pair_strip(cluster_list, horiz_center, half_width):
    'find closest pair in  the strip'
    newlist=[]
    # print '\n closest_pair_strip:---------------->',cluster_list,'<----------------\n'

    for item in cluster_list:
        if abs(item.horiz_center()-horiz_center)<=half_width:
            newlist.append(item)
    newlist.sort(key = lambda newlist: newlist.vert_center())
    # print '|||||....',newlist,horiz_center,half_width
    length_n=len(newlist)
    closest_pair=[float('inf'),-1,-1]
    if length_n-2==0:
        point1,point2=newlist[0],newlist[1]
        dist=point1.distance(point2)
        closest_pair[0]=dist

        closest_pair[1]=cluster_list.index(point1)

        closest_pair[2]=cluster_list.index(point2)
        if closest_pair[1]>closest_pair[2]:

            closest_pair[1],closest_pair[2]=closest_pair[2],closest_pair[1]
        return tuple(closest_pair)
    else:

        for u_item in xrange(length_n+1-2):
            limit=min(u_item+4,length_n+1-1)
            for v_item in xrange(u_item+1,limit):

                dist=newlist[u_item].distance(newlist[v_item])
                # print u_item,v_item,dist,range(u_item+1,limit)
                if closest_pair[0]>dist:
                    closest_pair[0]=dist

                    closest_pair[1]=cluster_list.index(newlist[u_item])
                    closest_pair[2]=cluster_list.index(newlist[v_item])
        if closest_pair[1]>closest_pair[2]:

            closest_pair[1],closest_pair[2]=closest_pair[2],closest_pair[1]
        # print '--',closest_pair
    return tuple(closest_pair)
def hierarchical_clustering(cluster_list, num_clusters):
    'find K clusters in points or data by hierarchical_clustering'
    k_cluster=len(cluster_list)

    while k_cluster>num_clusters:
        close_pair=fast_closest_pair(cluster_list)
        point1,pomit2=cluster_list[close_pair[1]],cluster_list[close_pair[2]]
        print point1,point1
        new_point=point1.merge_clusters(pomit2)
        cluster_list.append(new_point)
        cluster_list.remove(point1)
        cluster_list.remove(pomit2)
        k_cluster=len(cluster_list)
    return cluster_list
def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    'find K clusters in some points or data by hierarchical_clustering'
    total_number=len(cluster_list)
    init_centroids=cluster_list[:num_clusters]
    for n_time in range(num_iterations):
        emtpys=[]
        for m in range(num_clusters):
            emtpys.append([])
        for point in cluster_list:

                dist=map(lambda x:x.distance(point),init_centroids)
                dist.remove(0.0)
                dist.remove(0)
                ind=dist.index(min(dist))#第几个centroids最近
                emtpys[ind].extend([point,init_centroids[ind]])
        #updata cent
        newcent=list(init_centroids)
        for cent_index in range(num_clusters):
            mass=set(init_centroids[cent_index])
            point=reduce(lambda point1,ponit2:point1.merge_clusters(ponit2),mass)
            newcent[cent_index]=point
        init_centroids=list(newcent)
    return