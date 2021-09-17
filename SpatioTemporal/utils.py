import random
import matplotlib.pyplot as plt
from .models import *
from math import ceil

# Generate simulated samples in Saronic extend.
MIN_X_RANDOM = 23
MAX_X_RANDOM = 24
MIN_Y_RANDOM = 37
MAX_Y_RANDOM = 38
VARIANCE = 0.06

def _generate_traj(N : int, n_clusters : int, time_range : int, noise : bool = False) -> List[Trajectory]:
    """Generates simulated trajectories in the extend of Saronic.

    Parameters
    ----------
    N : int
        Number of trajectories to generate.
    n_clusters : int
        Number of clusters to simulate.
    time_range : int
        The maximum timestamp of the generated trajectory points, [0, time_range]
    noise : bool, default=False
        Generate noise if True.

    Returns
    -------
    trajs : list
        The generated trajectories.
    """
    if not noise: 
        # Create 'n_clusters' points as seed for the simulated clusters.
        X = [random.uniform(MIN_X_RANDOM,MAX_X_RANDOM) for _ in range(n_clusters)] 
        Y = [random.uniform(MIN_Y_RANDOM,MAX_Y_RANDOM) for _ in range(n_clusters)]

    trajs = []
    for _ in range(N):
        while True:
            # random_time = random.randint(0,3)
            # if random_time == 0:
            #     ts = time_range - random.randint(0, time_range)
            # else:
            #     ts = time_range - int(random.gauss(time_range/2, time_range/100))

            # Random generation of start time and end time of trajectory.
            ts = time_range - random.randint(0, time_range)
            te = ts + random.randint(0, time_range - ts)
            
            if te-ts != 0:
                break

        if not noise:
            # Generate start and end trajectory point "near" the random cluster seed.
            i = random.randint(0, n_clusters - 1)
            xs, xe = (random.gauss(X[i], VARIANCE) for _ in range(2))
            ys, ye = (random.gauss(Y[i], VARIANCE) for _ in range(2))
        else:
            # In case of noise generation, generate random point in the extend.
            x = random.uniform(MIN_X_RANDOM + VARIANCE, MAX_X_RANDOM - VARIANCE)
            y = random.uniform(MIN_Y_RANDOM + VARIANCE, MAX_Y_RANDOM - VARIANCE)
            xs, xe = (random.gauss(x, VARIANCE) for _ in range(2))
            ys, ye = (random.gauss(y, VARIANCE) for _ in range(2))

        tps = TPoint(xs, ys, ts)
        tpe = TPoint(xe, ye, te) 
        trajs.append(Trajectory(tps, tpe))

    return trajs


def generate_samples(N : int = 1000, p : float = 0.5, n_clusters : int = 3, time_range : int = 100, random_state : int = None) -> List[Trajectory]:
    """Generates 'N' trajectories where (1-p)*100 % is simulated in 'n_clusters' and p*100% is noise,
    having random start timestamp < end timestamp in [0, 'time_range'].

    Parameters
    ----------
    N : int, default=1000
        Number of trajectories to generate.
    p : float, default=0.3 
        Percentage of trajectories to generate as noise.
    n_clusters : int, default=3
        Number of clusters to simulate.
    time_range : int, default=100
        The maximum timestamp of the generated trajectory points, [0, time_range]
    random_state : int, default=None
        Random state to seed.

    Returns
    -------
    trajs : list
        List of clustered and noise simulated trajectories.
    """
    if random_state is not None:
        random.seed(random_state)

    # Generate clusters
    trajs = _generate_traj(int(N*(1-p)), n_clusters, time_range, noise=False)

    # Generate noise
    trajs_noise = _generate_traj(int(N*p), None, time_range, noise=True)

    return trajs + trajs_noise


def display(trajs : List[Trajectory] = None, 
            start_end : bool = False,
            links : bool = False, 
            mid : bool = False, 
            R : List[List[Cluster]] = None, 
            N : List[List[Cluster]] = None, 
            w : int = None, 
            V : List[List[GetisCluster]] = None, 
            ais : bool = False):
    """Display data.

    Parameters
    ----------
    trajs : list, default=None
        List of Trajectory objects to display.
    start_end : bool, default=False
        Display start and end points.
    links : bool, default=False
        Display the link between start point and end point. (requires 'trajs')
    mid : bool, default=False
        Display the midpoint of start point and end point. (requires 'trajs')
    R : list, default=None
        The R set to display. (requires N and w)
    N : list, default=None
        The N set to display. (requires R and w)
    w : int, default=None
        The time window of R and N. (requires R and N or GC)
    V : list
        A list of lists of GetisCluster objects having calculated Gi* and updated spot characterization ('Hot'/'Cold'/None) (requires w)
    ais : bool, default=False
        Trajs are ais data or not, specifies the extend.    
    """
    if trajs is not None:

        Xs = [tr.tps.x for tr in trajs]
        Ys = [tr.tps.y for tr in trajs]

        Xe = [tr.tpe.x for tr in trajs]
        Ye = [tr.tpe.y for tr in trajs]

        Xm = [tr.m.x for tr in trajs]
        Ym = [tr.m.y for tr in trajs]
    
        fig, ax = plt.subplots(1, 2, squeeze=False, figsize=(25,8))
        if start_end:
            ax[0][0].scatter(Xs, Ys, c='r', alpha=0.5, label='start')
            ax[0][0].scatter(Xe, Ye, c='b', alpha=0.5, label='end')
        if mid:
            ax[0][0].scatter(Xm, Ym, c='g', alpha=0.5, label='middle')
        if links:
            ax[0][0].plot([Xs, Xe], [Ys, Ye], c='c', alpha=0.5, label='link')
        ax[0][0].set_title("Data")
        ax[0][0].legend()
        
        durations = [tr.tpe.t - tr.tps.t for tr in trajs]
        ax[0][1].hist(durations, bins=40)
        ax[0][1].set_title("Duration")
        
        plt.show()
    elif R is not None and N is not None and w is not None:
        n = len(R)
        rows = ceil(n/3)
        fig, ax = plt.subplots(nrows=rows, ncols=3, figsize=(15,rows*5), squeeze=False)

        for i in range(rows):
            for j in range(3):
                index = i*3 + j
                if index >= n:
                    continue
                x = [p.x for cluster in R[index] for p in cluster.cpoints]
                y = [p.y for cluster in R[index] for p in cluster.cpoints]
                c = [cluster.c for cluster in R[index] for p in cluster.cpoints]
                clusters = ax[i][j].scatter(x, y, c=c, cmap='Set1', alpha=0.3)

                x = [p.x for cluster in N[index] for p in cluster.cpoints]
                y = [p.y for cluster in N[index] for p in cluster.cpoints]
                c = [i for i, cluster in enumerate(N[index]) for _ in cluster.cpoints]
                noise = ax[i][j].scatter(x, y, c=c, s=80, edgecolor='black', marker='*', cmap='tab10', label='noise', alpha=0.5)

                if ais:
                    extent = [23.03014957562005, 23.85354971429277, 37.7, 38.06] # Extend for ais2018db
                else:
                    extent = [MIN_X_RANDOM-2*VARIANCE, MAX_X_RANDOM+2*VARIANCE, MIN_Y_RANDOM-2*VARIANCE, MAX_Y_RANDOM+2*VARIANCE]
                ax[i][j].axis(xmin=extent[0], xmax=extent[1], ymin=extent[2], ymax=extent[3])

                ax[i][j].legend([clusters, noise], ['Clusters', 'Noise'])

                count_clusters = []
                for cluster in R[index]:
                    if cluster.c not in count_clusters:
                        count_clusters.append(cluster.c)
                count_clusters = len(count_clusters)
                count_noise = len(N[index])
                ax[i][j].set_title(f"Time slices: {index} -  {index+w-1}\n# clusters: {count_clusters}, # noise: {count_noise}")
        fig.tight_layout(pad=0.3)    
        plt.show()
    elif V is not None and w is not None:
        n = len(V)
        rows = ceil(n/3)
        fig, ax = plt.subplots(nrows=rows, ncols=3, figsize=(15,rows*5), squeeze=False)

        for i in range(rows):
            for j in range(3):
                index = i*3 + j
                if index >= n:
                    continue
                hx = []
                hy = []
                cx = []
                cy = []
                nx = []
                ny = []
                count_hot = 0
                count_cold = 0
                count_none = 0
                for gc in V[index]:
                    if gc.spot == 'Hot':
                        hx += [cp.x for cp in gc.cpoints]
                        hy += [cp.y for cp in gc.cpoints]
                        count_hot += 1
                    elif gc.spot == 'Cold':
                        cx += [cp.x for cp in gc.cpoints]
                        cy += [cp.y for cp in gc.cpoints]
                        count_cold += 1
                    else:
                        nx += [cp.x for cp in gc.cpoints]
                        ny += [cp.y for cp in gc.cpoints]
                        count_none += 1
                ax[i][j].scatter(nx, ny, c='k', alpha=0.1)
                ax[i][j].scatter(hx, hy, c='r', alpha=0.1)
                ax[i][j].scatter(cx, cy, c='b', alpha=0.1)

                if ais:
                    extent = [23.03014957562005, 23.85354971429277, 37.7, 38.06] # Extend for ais2018db
                else:
                    extent = [MIN_X_RANDOM-2*VARIANCE, MAX_X_RANDOM+2*VARIANCE, MIN_Y_RANDOM-2*VARIANCE, MAX_Y_RANDOM+2*VARIANCE]

                ax[i][j].axis(xmin=extent[0], xmax=extent[1], ymin=extent[2], ymax=extent[3])

                ax[i][j].set_title(f"Time slices: {index} - {index+w-1}\n# hot: {count_hot}, # cold {count_cold}, # none {count_none}")

                # Uncomment this and comment hx=[] until .set_title() to get the histograms of xs
                # values = [gc.x for gc in V[index]]
                # ax[i][j].hist(values, bins=20)
                # avg_values = sum(values)/len(values)
                # ax[i][j].plot([avg_values, avg_values], [0, 30])
        fig.tight_layout(pad=0.3)
        plt.show()

