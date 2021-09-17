from typing import List

class Point:
    """Point class.

    Parameters
    ----------
    x : float
        x coordinate.
    y : float
        y coordinate.
    """
    def __init__(self, x : float, y : float):
        self.x = x
        self.y = y

class TPoint(Point):
    """TPoint class, i.e. Trajectory Point.

    Parameters
    ----------
    x : float
        x coordinate.
    y : float
        y coordinate.
    t : int
        Timestamp.
    """
    def __init__(self, x : float, y : float, t : int):
        super().__init__(x, y)
        self.t = t


class Trajectory:
    """Trajectory class, expresses the start point and end point of the 'skein'.

    Parameters
    ----------
    tps : TPoint
        The starting trajectory point.
    tpe : TPoint
        The ending trajectory point.
    m : Point
        The centroid of tps and tpe
    """
    def __init__(self, tps, tpe, m = None):
        self.tps = tps
        self.tpe = tpe
        self.m = Point((tpe.x + tps.x)/2, (tpe.y + tps.y)/2) if m is None else m 

class TPartition:
    """TPartition class, the timeslice with its corresponding points.

    Parameters
    ----------
    ts : int
        The lower time limit of partition.
    te : int 
        The upper time limit of partition.
    tpoints : list
        A list of TPoints in partition.
    """
    def __init__(self, ts : int, te : int, tpoints : List[TPoint]):
        self.ts = ts
        self.te = te
        self.tpoints = tpoints


class ClusterPoint(TPoint):
    """ClusterPoint class, i.e TPoint in a cluster or just noise (c = -1)

    Parameters
    ----------
    x : float
        x coordinate.
    y : float
        y coordinate.
    t : int
        Timestamp.
    c : int
        Cluster label.
    """
    def __init__(self, x : float, y : float, t : int, c : int):
        super().__init__(x, y, t)
        self.c = c

class Cluster:
    """Cluster class, a cluster of cluster points.

    Parameters
    ----------
    id : int
        A unique identifier.
    cpoints : list
        A list of ClusterPoint objects.
    t : int
        Timestamp
    c : int, default=None
        Cluster label
    """
    def __init__(self, id : int, cpoints : List[ClusterPoint], t : int, c : int = None):
        self.id = id
        self.cpoints = cpoints
        mx, my = sum([cp.x for cp in cpoints]) / len(cpoints), sum([cp.y for cp in cpoints]) / len(cpoints)
        self.m = Point(mx, my)
        self.t = t
        self.c = c

class GetisCluster(Cluster):
    """GetisCluster class, the cluster of points to evaluate its Gi* and statisticaly characterize as 'Hot', 'Cold' or None.

    Parameters
    ----------
    id : int
        A unique identifier.
    cpoints : list
        A list of ClusterPoint objects.
    t : int
        Timestamp

    Attributes
    ----------
    m : Point
        The centroid to use in distance calculations.
    x : int
        The feature of interest, that is the number of ClusterPoints of the GetisCluster.
    gi : float
        The Gi* statistic 
    significant : bool
        A statistical significance flag, i.e. it is hot spot or cold spot.
    spot : str
        The characterization of the cluster, 'Hot', 'Cold' or None.    
    """
    def __init__(self, id : int, cpoints : List[ClusterPoint], t : int):
        super().__init__(id = id, cpoints = cpoints, t = t)
        self.x = len(cpoints)
        self.gi = None
        self.significant = False
        self.spot = None