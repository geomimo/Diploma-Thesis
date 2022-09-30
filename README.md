# Detecting statistically significant spatial clusters in vessel data

## Diploma Thesis

Author: *George Mimoglou*  
Supervisor: *Yannis Theodoridis, Professor*  
Location: *University of Piraeus - 2021*

## Run Tests

- Create a conda enviroment using the `environment.yml`
- Execute all cells in notebooks `test_ais.ipynb` and `test_dummy.ipynb`

## Overview

In the scope of the dissertation, two methods were implemented and tested on simulated and real data. Both detect spatial clusters (anchorages), but each one characterizes them either as *long/short-term* based on their temporal persistence `Spatio-Temporal Clustering (STC)`, or as *hot/cold spot* based on their spatio-temporal autocorrelation `Spatio-Temporal Hot Spots (STHS)`. These can answer questions like *"how do anchorages influence each other in the spatio-temporal dimensions?"*. Both methods got inspired by the grid-based method `Trajectory Hot Spots`[^1], where cells in a spatio-temporal partition are clustered and statistically characterized based on the total duration of vessel existence.

Most statistical spatial analysis methods of modern tools use *Complete Spatial Randomness (CSR)* as the Null Hypothesis, i.e., spatial data does not present any underline structure. Although they excel when applied to dense data, e.g., states in a country, they produce misleading results on sparse grids of maritime data. Furthermore, the results of grid-based methods depend significantly on size parameters and are subject to the Modifiable Areal Unit Problem[^2].

The proposed methods overcome these issues by using a series of clustering in a rolling window manner and partitioning only the time dimension with domain-interpretable parameters.

The methods were tested using simulated and real AIS data; see `Run tests`.

### Preprocessing

Even though both methods apply to various domains, anchorage cases require specific data preprocessing. Vessel data are [Automated Identification System (AIS)](https://www.imo.org/en/OurWork/Safety/Pages/AIS.aspx) signals, providing position, speed, rate of turn, heading, and vessel identity, among others, so it is needed to transform them into an appropriate form to extract anchorages.

In this implementation, the points indicating that the vessel enters and exits a possible anchorage are extracted using simple heuristics, for example, very low speed and rapid change of bearing. Then, their midpoint is calculated, indicating a stationary point in a possible anchorage. Lastly, we partition time and assign each midpoint to the timeslices it exists. The models are fitted to a list of midpoints.

AIS data collected from the Argosaronic Gulf were used for testing. The noise was cleansed using mask operations in QGIS and spatial queries in PostgreSQL.

### Spatio-Temporal Clustering (STC)

This model breaks down the meaning of spatio-temporal autocorrelation using logical reasoning and tries to approximate it without using statistics. Intuitively, the sparsity of the data and the stationarity of anchorages allow us to focus on temporal characteristics. In other words, we expect that anchorages do not influence each other in the spatial domain. However, each instance can influence itself in the future, thus existing for either long or short durations.

In the first phase, STC applies a clustering algorithm, `g_clustering`, on the midpoints (see `Preprocessing`) for each timeslice to detect the anchorages for a particular period (timeslice). The noise should be ignored. In the second phase, the centroid for each cluster is calculated. Then, a second clustering-with-noise algorithm, `v_clustering`, is applied to the centroids of $w$ timeslices projected onto a common plane in a rolling window manner. The detected clusters are the long-term anchorages, and the noise is the short-term anchorages.

The first clustering approximates the spatial autocorrelation notion by detecting structure, i.e., anchorages, while the second clustering integrates time and encapsulates the association of an anchorage with itself through time.

The choice of the clustering algorithms is free; however, they are associated with interpreting the results. `g_clustering` is used to detect anchorages, so KMeans, which requires prior knowledge of clusters, is not an appropriate choice. In contrast, with DBSCAN, one can determine the minimum number of vessels `n_samples` to exist in a proximity `eps` to detect an anchorage. At `v_clustering`, DBSCAN parameters can be used to determine the minimum duration or instances of an anchorage  `n_samples`, out of $w$ timeslices to characterize it as long-term or short-term otherwise.

### Spatio-Temporal Hot Spots (STHS)

STHS follows the same rationale, but the second clustering, `v_clustering` is replaced with hot spot analysis using *Getis-Ord* $G^*_i$[^3]. It differentiates itself from STC since it considers the cardinality of the clusters. STC aggregates this information into a single centroid without paying attention to the evolution of cardinality. For example, a long-term anchorage can be a hot spot if there is high cardinality for a couple of consecutive timeslices. There is spatial autocorrelation due to the stationarity of the cluster, and temporal autocorrelation due to the consistently high cardinality. Also, after some timeslices, the same anchorage can be a cold spot when the number of vessels declines.

Once `g_clustering` finishes, each cluster midpoint and its cardinality are calculated. Then, using the same rolling window method, the $G^*_i$ is calculated considering cardinality as the characteristic attribute. *Getis-Ord* $G^*_i$ requires the spatial weight between two clusters. For the spatio-temporal dimensions, we have defined the *Inverse Distance Weighting (IDW)* as:

$$w_{K_{xa}, K_{yb}} = {m \over e^{\lparen z+1 \rparenφ\lparen K_{xa}, K_{yb}\rparen

where $K_{xa}, K_{yb}$ anchorages from timeslices $x$ and $y$, $z=\lvert x-y \rvert$ the distance of the timeslices and $φ$ a spatial distane function. $m$ is related to scale.

The clusters are characterized as `Hot`, `Cold` spots, or `None` based on the *G^\*^i* value and the sign of the *z-scores*.

## Keywords

Maritime Surveillance, Clustering, Statistically Significant Regions, Anchorages, AIS, Spatio-temporal Dimensions

[^1]: P. Nikitopoulos, A. Paraskevopoulos, C. Doulkeridis, N. Pelekis and Y. Theodoridis, "Hot Spot Analysis over Big Trajectory Data," 2018 IEEE International Conference on Big Data (Big Data), 2018, pp. 761-770, doi: 10.1109/BigData.2018.8622376.

[^2]: Openshaw, Stan (1984). The modifiable areal unit problem. Norwick: Geo Books. ISBN 0860941345. OCLC 1205248

[^3]: Griffith, Daniel. (1992). What is Spatial Autocorrelation? Reflections on the Past 25 Years of Spatial Statistics. Espace géographique. 21. 265-280. 10.3406/spgeo.1992.3091.
