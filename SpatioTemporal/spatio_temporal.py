from matplotlib.pyplot import vlines
from sklearn.metrics.pairwise import haversine_distances
from .models import * 
from sklearn.cluster import DBSCAN
from collections import OrderedDict
import numpy as np
from math import sqrt, exp, radians
from typing import Dict, List, Literal, Tuple
from sklearn.base import BaseEstimator


class ST:
    """Base class of STC and STHS.

        Parameters
        ----------
        g_clustering : BaseEstimator, default=DBSCAN
            A class that perfroms clustering, having the .fit() method. It is used for the 
            clustering of each Gj.
        g_clustering_args : dict, default={'eps': 1, 'min_samples': 2}
            Dict-like args for the g_clustering.
        interval : int, default=10
            The interval of timeslices to create (T time partition).
        w : int, default=3
            The number of time partitions (time window) to aggregate and perfom v_clustering. 
    """
    def __init__(self, g_clustering : BaseEstimator = DBSCAN, 
                       g_clustering_args : Dict[str, object] = {'eps': 1, 'min_samples': 2}, 
                       interval : int = 10,
                       w : int = 3):
        self.g_clustering = g_clustering
        self.g_clustering_args = g_clustering_args
        self.interval = interval
        self.w = w


    def _create_G(self, trajs : List[Trajectory]) -> List[TPartition]:
        """Creates the G set, i.e. the set of sets of points for each time partition.
        G = {G_1, G_2, ..., G_n}
        G_1 = {g_11, g_12, ..., g_1m}

        Paramters
        ---------
        trajs : list
            A list of Trajectory objects
        
        Returns
        -------
        G : list
            A list of TPartition objects, i.e. the Gj, starting time, end time for that partition.
        """
        G = []
        time_range = max(tr.tpe.t for tr in trajs)
        for t in range(0, time_range, self.interval):
            t1 = t
            t2 = t + self.interval
            Gj = []
            for tr in trajs:
                if tr.tps.t > t2 or tr.tpe.t < t1:
                    continue
                m_new = TPoint(tr.m.x, tr.m.y, t1)

                Gj.append(m_new)
            G.append(TPartition(t1, t2, Gj))

        self.G = G
        return G


    def _g_clustering(self, G : List[TPartition]) -> List[List[Cluster]]:
        """Perfoms clustering for every Gj set, i.e. m points of each time partition.

        Parameters
        ----------
        G : list 
            A list of TPartition objects, i.e. the Gj, starting time, end time for that partition.

        Returns
        -------
        K : list
            A list of lists of Cluster objects, i.e. the clusters of each time partition.
        """
        K = []
        id = 0
        for Gj in G:
            data = np.array([[m.x, m.y] for m in Gj.tpoints])

            if len(data) == 0:
                K.append([])
                continue

            model = self.g_clustering(**self.g_clustering_args)
            model.fit(data)            

            clusters = {}
            for i in range(len(data)):
                if model.labels_[i] == -1: # ignore noise
                    continue

                if model.labels_[i] not in clusters:
                    clusters[model.labels_[i]] = []
                clusters[model.labels_[i]].append(data[i])

            Kj = []
            for k, v in clusters.items():
                cps = []
                for coords in v:
                   cps.append(ClusterPoint(coords[0], coords[1], Gj.ts, k))
                Kj.append(Cluster(id = id, cpoints = cps, t = Gj.ts, c = k))
                id += 1
            K.append(Kj)
        return K


class STC(ST):
    """Performs Spatio-Temporal Clustring from list of Trejectory objects.

    Parameters
    ----------
    g_clustering : BaseEstimator, default=DBSCAN
        A class that perfroms clustering, having the .fit() method. It is used for the 
        clustering of each Gj.
    g_clustering_args : dict, default={'eps': 1, 'min_samples': 2}
        Dict-like args for the g_clustering.
    v_clustering : BaseEstimator, default=DBSCAN
        A class that perfroms clustering, having the .fit() method. It is used fot the 
        clustring of each Vj.
    v_clustering_args : dict, default={'eps': 1, 'min_samples': 2}
        Dict-like args for the v_clustering.
    interval : int, default=10
        The interval of timeslices to create (T time partition).
    w : int, default=3
        The number of time partitions (time window) to aggregate and perfom v_clustering. 
    """
    def __init__(self, g_clustering : BaseEstimator = DBSCAN, g_clustering_args : Dict[str, object] = {'eps': 1, 'min_samples': 2},
                       v_clustering : BaseEstimator = DBSCAN, v_clustering_args : Dict[str, object] = {'eps': 1, 'min_samples': 2},
                       interval : int = 10, w : int = 3):
        super().__init__(g_clustering=g_clustering, g_clustering_args=g_clustering_args, interval=interval, w = w)
        self.v_clustering = v_clustering
        self.v_clustering_args = v_clustering_args
        
    

    def _v_clustering(self, K : List[List[Cluster]]) -> Tuple[List[List[Cluster]], List[List[Cluster]]]:
        """Perfoms clustering for every Vj set, i.e. union of w Kjs
        
        Parameters
        ----------
        K : list
            A list of lists of Cluster objects, i.e. the clusters of each time partition.

        Returns
        -------
        R : list
            A list of lists of Cluster objects with updated cluster label, i.e. the clusters of each Vj.
        N : list
            A list of lists of Cluster object with cluster label equal to -1 (noise), i.e. the noise of each Vj.
        """
        R = []
        N = []
        for i in range(len(K) - self.w + 1):

            Vj = [K[j] for j in range(i, i + self.w)]
            Vj = [item for sublist in Vj for item in sublist] # flatten list

            x1 = [cluster.m.x for cluster in Vj]
            x2 = [cluster.m.y for cluster in Vj]
            data = list(zip(x1,x2))
            if len(data) == 0:
                R.append([])
                N.append([])
                continue
            
            model = self.v_clustering(**self.v_clustering_args)
            model.fit(data)

            Rj = [Cluster(i, cpoints=cluster.cpoints, t=cluster.t, c=model.labels_[i]) for i, cluster in enumerate(Vj) if model.labels_[i] != -1]
            Nj = [Cluster(i, cpoints=cluster.cpoints, t=cluster.t, c=model.labels_[i]) for i, cluster in enumerate(Vj) if model.labels_[i] == -1]

            R.append(Rj)
            N.append(Nj)

        return R, N


    def fit(self, trajs : List[Trajectory]) -> Tuple[List[List[Cluster]], List[List[Cluster]]]:
        """Performs Spatio-Temporal Clustring from list of Trejectory objects.

        Paramters
        ---------
        trajs : list
            A list of Trajectory objects

        Returns
        -------
        R : list
            A list of lists of Cluster objects with updated cluster label, i.e. the clusters of each Vj.
        N : list
            A list of lists of Cluster object with cluster label equal to -1 (noise), i.e. the noise of each Vj.
        """
        G = self._create_G(trajs)
        K = self._g_clustering(G)
        R, N = self._v_clustering(K)
        return R, N 


class STHS(ST):
    """Performs Spatio-Temporal Hot Spot from list of Trejectory objects.

    Notes
    -----
    Haversine distance is used as metric for the weight function.
    m in weight function is set to 1.

    Parameters
    ----------
    g_clustering : BaseEstimator, default=DBSCAN
        A class that perfroms clustering, having the .fit() method. It is used for the 
        clustering of each Gj.
    g_clustering_args : dict, default={'eps': 1, 'min_samples': 2}
        Dict-like args for the g_clustering.
    interval : int, default=10
        The interval of timeslices to create (T time partition).
    w : int, default=3
        The number of time partitions (time window) to aggregate and calculate Getis-Ord Gi*. 
    """
    def _create_getis_clusters(self, K : List[List[Cluster]]) -> List[List[GetisCluster]]:
        """Creates a list of lists of GetisCluster objects from a list of lists of GetisCluster objects.

        Parameters
        ----------
        K : list
            A list of lists of Cluster objects, i.e. the clusters of each time partition.

        Returns
        -------
        GC : list
            A list of lists of GetisCluster objects, i.e the getis clusters if each time partition.
        """
    
        GC = []
        for Kj in K:
            GCj = []
            for cluster in Kj:
                GCj.append(GetisCluster(cluster.id, cluster.cpoints, cluster.c))
            GC.append(GCj)
        return GC


    def _calculate_distance(self, GC : List[List[GetisCluster]]) -> Dict[int, Dict[int, float]]:
        """Calculates a distance matrix of all GetisCluster objects.
        
        Notes
        -----
        Haverinse distance is calculated in this implementation.

        Parameters
        ----------
        GC : list
           A list of lists of GetisCluster objects, i.e the getis clusters if each time partition.

        Returns
        -------
        distances : dict
            A dict of dicts with clusters' id as keys and the distance as value
        """
        distances = {}
        for i in range(len(GC)):
            GCja = GC[i]
            for clustera in GCja:
                distances[clustera.id] = {}
                p1 = [[radians(clustera.m.y), radians(clustera.m.x)]] #(lat, lon)
                try:
                    for j in range(i, i + self.w):
                        GCjb = GC[j]
                        for clusterb in GCjb:
                            p2 = [[radians(clusterb.m.y), radians(clusterb.m.x)]]
                            distances[clustera.id][clusterb.id] = haversine_distances(p1,p2)[0][0]
                except IndexError:
                    continue

        return distances


    def _calculate_getis_ord(self, GC : List[List[GetisCluster]], a : float = 0.01) -> List[List[GetisCluster]]:
        """Calculates the Gi* statistic for every GetisCluster in Vj, for every Vj.
        
        Notes
        -----
        Haversine distance is used as metric for the weight function.
        m in weight function is set to 1. 

        Parameters
        ----------
        GC : list
           A list of lists of GetisCluster objects, i.e the getis clusters if each time partition.
        a : int, default=0.10
            The significance level for the hypothesis testing.
        
        Returns 
        -------
        V : list
            A list of lists of GetisCluster objects having calculated Gi* and updated spot characterization ('Hot'/'Cold'/None)
        """
        distances = self._calculate_distance(GC)

        def weight(z, gc1, gc2):
            try:
                return 1*exp(-(z+1)*distances[gc1.id][gc2.id])
            except KeyError:
                return 1*exp(-(z+1)*distances[gc2.id][gc1.id])
            

        V = []
        for i in range(len(GC) - self.w + 1):

            Vj = [gc for j in range(i, i + self.w) for gc in GC[j]]           
            n = len(Vj)
            if n == 0:
                V.append(Vj)
                continue
            x_bar = sum(gc.x for gc in Vj) / n
            S = sqrt(sum(gc.x ** 2 for gc in Vj) / n - x_bar ** 2) 

            for gc1 in Vj:
                swx = 0
                sw = 0
                sw2 = 0
                for gc2 in Vj:
                    z = abs(gc1.t - gc2.t)
                    w = weight(z, gc1, gc2)
                    swx += w * gc2.x  
                    sw += w
                    sw2 += w ** 2
                gc1.gi = (swx - x_bar * sw) / (S * sqrt((n * sw2 - sw ** 2) / (n - 1)))
                if a == 0.10:
                    critical = 1.65
                elif a == 0.05:
                    critical = 1.96
                elif a == 0.01:
                    critical = 2.58
                else:
                    raise ValueError("'a' must be 0.10, 0.05 or 0.01.")
                
                if gc1.gi <= -critical:
                    gc1.significant = True
                    gc1.spot = 'Cold'
                elif gc1.gi >= critical:
                    gc1.significant = True
                    gc1.spot = 'Hot'
            V.append(Vj)
            
        return V
                

    def fit(self, trajs : List[Trajectory], a : float = 0.05) -> List[List[GetisCluster]]:
        """Performs Spatio-Temporal Hot Spot from list of Trejectory objects.
        
        Parameters
        ----------
        trajs : list
            A list of Trajectory objects
        a : int, default=0.10
            The significance level for the hypothesis testing.

        Returns
        -------
        V : list
            A list of lists of GetisCluster objects having calculated Gi* and updated spot characterization ('Hot'/'Cold'/None)
        """
        G = self._create_G(trajs)
        K = self._g_clustering(G)
        GC = self._create_getis_clusters(K)
        V = self._calculate_getis_ord(GC, a)
        return V



