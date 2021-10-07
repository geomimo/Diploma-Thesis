from SpatioTemporal.utils import generate_samples, display
from SpatioTemporal.spatio_temporal import STC, STHS
from SpatioTemporal.models import TPoint, Trajectory, Point
from tabulate import tabulate
import pandas as pd
import geopandas as gpd
from preprocessing import preprocessing
from shapely import wkt
from typing import List, Dict

def _procedure(trajs : List[Trajectory], time_range : int, interval : int, w : int, 
               g_clustering_args : Dict[str, object] = {'eps':1, 'min_samples': 2}, 
               v_clustering_args : Dict[str, object] = {'eps':1, 'min_samples': 3}, 
               ais : bool = True, a : float = 0.05):
    """Executes STC and STHS and displays the results.
    """    
    stc = STC(interval=interval, w=w, g_clustering_args=g_clustering_args, v_clustering_args=v_clustering_args)
    
    R, N = stc.fit(trajs)
    display(trajs=trajs, start_end=False, links=False, mid=True)
    display(R=R, N=N, w=w, ais=ais)

    sths = STHS(interval=interval, w=w, g_clustering_args=g_clustering_args)
    V = sths.fit(trajs, a)
    display(V=V, w=w, ais=ais)


def _to_trajectories(df) -> List[Trajectory]:
    """Transforms dataframe of ms (ms.csv calculated in preprocessing.py) in Trajectory list.
    """
    trajs = []
    for i, row in df.iterrows():
        tps = TPoint(None, None, row['t_s'])
        tpe = TPoint(None, None, row['t_e'])
        m = Point(row.geom.coords[0][0], row.geom.coords[0][1])
        trajs.append(Trajectory(tps, tpe, m=m))

    return trajs


def _load_ms() -> gpd.GeoDataFrame:
    """Loads the ms.csv (calculated in preprocessing.py) in GeoDataFrame.
    """
    m_df = pd.read_csv('./ms.csv')
    m_df.geom = m_df.geom.apply(wkt.loads)
    m_df = gpd.GeoDataFrame(m_df, crs=4326, geometry='geom')
    m_df['t_s'] = pd.to_datetime(m_df['t_s'])
    m_df['t_e'] = pd.to_datetime(m_df['t_e'])
    return m_df


def _convert_timestamp(df) -> gpd.GeoDataFrame:
    """Converts timestamps to intergers in [0, max_timestamp - min_timestamp] range.
    """
    df['t_s'] = pd.to_datetime(df['t_s'])
    df['t_e'] = pd.to_datetime(df['t_e'])

    df['t_s'] = df['t_s'].view('int64') // 10 ** 9
    df['t_e'] = df['t_e'].view('int64') // 10 ** 9
    min_t = df['t_s'].min()
    df['t_s'] -= min_t
    df['t_e'] -= min_t

    return df


def real_test(i=-1, preprocess=False, interval=None, w=None, g_clustering_args=None, v_clustering_args=None, a=None):
    """Executes STC and STHS using ais data. If 'preprocess' == True then it computes the ms from raw data,
    else loads the precomputed ms from the ms.csv.

    Notes
    -----
    'preprocess' = True is slow, 15-30 min.  

    Parameters
    ----------
    i : int, default=-1
        Test to execute, available: 0, 1, 2, 3 
    preprocess : bool, default=False
        If True, executes preprocessing of AIS data. Requires 'dbcon.json' configuration.
    interval : int, default=None
        The interval of timeslices to create (T time partition).
    w : int, default=None
        The number of time partitions (time window) to aggregate.
    g_clustering_args : dict, default=None
        Dict-like args for the g_clustering. 
    v_clustering_args : dict, default=None
        Dict-like args for the v_clustering.
    a : float, default=None
        The significance level for the hypothesis testing.
    """
    trajs = []
    if preprocess:
        _, _, m_df = preprocessing()
    else:
        m_df = _load_ms()
    
    m_df = _convert_timestamp(m_df)

    time_range = m_df['t_e'].max()
    if i == 0:
        interval=60*60*24*2 # 2 days
        w = 15 # 2 days * 15 = 30 days
        g_clustering_args = {"eps": 0.010, "min_samples": 4}
        v_clustering_args = {"eps": 0.010, "min_samples": 5}
        a = 0.01
    elif i == 1:
        interval=60*60*24*2 # 2 days
        w = 15 # 2 days * 15 = 30 days
        g_clustering_args = {"eps": 0.010, "min_samples": 2}
        v_clustering_args = {"eps": 0.010, "min_samples": 10}
        a = 0.01
    elif i == 2:
        interval=60*60*24*7 
        w = 4
        g_clustering_args = {"eps": 0.013, "min_samples": 4}
        v_clustering_args = {"eps": 0.020, "min_samples": 4}
        a = 0.01
    elif i == 3:
        interval=60*60*24*31 
        w = 12
        g_clustering_args = {"eps": 0.020, "min_samples": 10}
        v_clustering_args = {"eps": 0.030, "min_samples": 12}
        a = 0.01
        

    trajs = _to_trajectories(m_df)

    data = [[len(trajs), time_range, interval, w, g_clustering_args, v_clustering_args, a]]
    headers = ["N", "time_range (s)", "interval (s)", "w (time partitions)", "g_clustering_args", "v_clustering_args", "a"]
    print(tabulate(data, headers=headers))
    _procedure(trajs, time_range, interval, w, g_clustering_args, v_clustering_args, ais=True) 


def dummy_test(i=-1, N=None, n_clusters=None, time_range=None, interval=None, w=None, random_state=None, g_clustering_args=None, v_clustering_args=None, a=None):
    """Executes STC and STGS using simulated data, generated by SpatioTemporal.utils.generate_samples.
    """
    
    if i == 0:
        N = 5000 
        n_clusters = 5 
        time_range = 300 
        interval = 15
        w = 4
        g_clustering_args = {"eps": 0.025, "min_samples": 10}
        v_clustering_args = {"eps": 0.040, "min_samples": 3}
        a = 0.01
        random_state = 21
    elif i == 1: 
        N = 1000
        n_clusters = 3
        time_range = 300
        interval = 5 
        w = 10
        g_clustering_args = {"eps": 0.010, "min_samples": 3}
        v_clustering_args = {"eps": 0.050, "min_samples": 6}
        a = 0.01
        random_state = 0
    elif i == 2:
        N = 10000
        n_clusters = 2
        time_range = 600
        interval = 100
        w = 3 
        g_clustering_args = {"eps": 0.015, "min_samples": 6}
        v_clustering_args = {"eps": 0.010, "min_samples": 2}
        a = 0.01
        random_state = 1
    elif i == 3:
        N = 5000
        n_clusters = 5
        time_range = 100
        interval = 1
        w = 30 
        g_clustering_args = {"eps": 0.020, "min_samples": 10}
        v_clustering_args = {"eps": 0.040, "min_samples": 25}
        a = 0.01
        random_state = 3

    data = [[N, n_clusters, time_range, interval, w, g_clustering_args, v_clustering_args, a]]
    headers = ["N", "n_clusters", "time_range (days)", "interval (days)", "w (time partitions)", "g_clustering_args", "v_clustering_args", "a"]
    print(tabulate(data, headers=headers))
    trajs = generate_samples(N=N, p=0.4, n_clusters=n_clusters, time_range=time_range, random_state=random_state)
    _procedure(trajs, time_range, interval, w, g_clustering_args, v_clustering_args, ais=False, a=a) 
