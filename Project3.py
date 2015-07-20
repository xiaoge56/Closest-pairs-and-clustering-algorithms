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
    closest_pair=[float('inf'),-1,-1]
    for u_item in xrange(len(cluster_list)):

        for v_item in xrange(u_item+1,len(cluster_list)):
            distance=cluster_list[u_item].distance(cluster_list[v_item])
            if closest_pair[0]>distance:
                closest_pair[0]=distance
                closest_pair[1]=u_item
                closest_pair[2]=v_item
    return tuple(closest_pair)

def fast_closest_pair(cluster_list):
    'find closest pair in cluster_list ,quick than slow_closest_pair()'
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
            print mid1,cluster_list,right_pair[1]+mid1
            closest_pair[1]=cluster_list.index(templist[right_pair[1]+mid1])
            closest_pair[2]=cluster_list.index(templist[right_pair[2]+mid1])
        mid2=1/2*(templist[mid1-1].horiz_center()+templist[mid1].horiz_center())
        other_point=closest_pair_strip(cluster_list,mid2,closest_pair[0])
        if closest_pair[0]>other_point[0]:
            closest_pair[0]=other_point[0]
            closest_pair[1]=cluster_list.index(other_point[1])
            closest_pair[2]=cluster_list.index(other_point[2])
    return tuple(closest_pair)
def closest_pair_strip(cluster_list, horiz_center, half_width):
    'find closest pair in  the strip'
    newlist=[]


    for item in cluster_list:
        if abs(item.horiz_center()-horiz_center)<=half_width:
            newlist.append(item)
    newlist.sort(key = lambda newlist: newlist.horiz_center())
    length_n=len(newlist)
    closest_pair=[float('inf'),-1,-1]
    if length_n-2==0:
        point1,point2=newlist[0],newlist[1]
        dist=point1.distance(point2)
        closest_pair[0]=dist

        closest_pair[1]=cluster_list.index(point1)

        closest_pair[2]=cluster_list.index(point2)

        return tuple(closest_pair)
    else:
        print newlist
        for u_item in xrange(length_n+1-2):
            limit=min(u_item+3,length_n+1-1)
            for v_item in xrange(u_item+1,limit):
                print u_item,v_item,u_item+1,limit
                dist=newlist[u_item].distance(newlist[v_item])
                if closest_pair[0]>dist:
                    closest_pair[0]=dist
                    closest_pair[1]=cluster_list.index(newlist[u_item])
                    closest_pair[2]=cluster_list.index(newlist[v_item])
    return tuple(closest_pair)
# m=closest_pair_strip([Cluster(set([]), 0, 0, 1, 0), Cluster(set([]), 1, 0, 1, 0), Cluster(set([]), 2, 0, 1, 0), Cluster(set([]), 3, 0, 1, 0)], 1.5, 1.0)
# m=closest_pair_strip([Cluster(set([]), 0.23, 0.94, 1, 0), Cluster(set([]), 0.65, 0.08, 1, 0),Cluster(set([]), 0.66, 0.43, 1, 0), Cluster(set([]), 0.91, 0.6, 1, 0), Cluster(set([]), 0.94, 0.9, 1, 0)], 0.65500000000000003, 0.30149599999999999)
# print 0.65500000000000003+0.30149599999999999,0.65500000000000003-0.30149599999999999
# print math.sqrt((0.91-0.66)**2+(0.6-0.43)**2),math.sqrt((0.03)**2+(0.3)**2)
# m=closest_pair_strip([Cluster(set([]), -1.0, 0.0, 1, 0),Cluster(set([]), -0.99, -10.0, 1, 0),Cluster(set([]), -0.98, -20.0, 1, 0), Cluster(set([]), 0.98, 20.0, 1, 0), Cluster(set([]), 0.99, 10.0, 1, 0), Cluster(set([]), 1.0, 0.0, 1, 0)], 0.0, 10.000005)
print fast_closest_pair([Cluster(set([]), 1.0, 0.0, 1, 0), Cluster(set([]), 4.0, 0.0, 1, 0), Cluster(set([]), 5.0, 0.0, 1, 0), Cluster(set([]), 7.0, 0.0, 1, 0)])  #expected one of the tuples in set([(1.0, 1, 2), (1.0, 0, 1), (1.0, 2, 3)])