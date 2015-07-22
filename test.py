# coding=utf-8
def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    'find K clusters in some points or data by hierarchical_clustering'
    total_n=num_clusters
    new_cluster_list=[c.copy() for c in cluster_list]
    new_cluster_list.sort(key = lambda new_cluster_list: new_cluster_list.total_population())

    centroids=new_cluster_list[-total_n:]
    print centroids

    fanhui=[]
    while num_iterations>0:
        num_iterations-=1
        retset=[]
        for m_num in range(num_clusters):
            retset.append([])
        for point in new_cluster_list:#开始对每个点进行染色
            dist=[]#每轮迭代中，点最靠近哪个类心
            poi=point.copy()
            for centroid in centroids:

                dist.append(point.distance(centroid))
            index=dist.index(min(dist))
            retset[index].append(poi)

        for index in range(num_clusters):
            for point in retset:
                if len(point)==1:
                    new_cent=point
                else:
                    # massset=[c.copy() for c in point]
                    massset=[c.copy() for c in point]
                    # massset=emtpys[cent_index]
                    new_cent=reduce(lambda x,y:x.merge_clusters(y),massset)

                    #
                    # list_del=[c.copy() for c in massset[1:]]
                    # for del_point in list_del:
                    #     pass
            centroids[index]=new_cent
        fanhui=retset
    real_list=[]
    for m in fanhui:
        if len(m)==1:
            real_list.append(m[0])
        else:
            a=m[0]
            for x in m[1:]:
                a.merge_clusters(x)
            real_list.append(a)




    return real_list