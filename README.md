# Detecting statistically significant spatial clusters in vessel data

## Diploma Thesis

Author: *George Mimoglou*<br>
Supervisor: *Yannis Theodoridis, Professor*<br>
Location: *University of Piraeus - 2021*

## Abstract
Surveillance of marine areas is major pillar of the modern maritime evolution, in both European and International level. This is due to the complex nature of the domain which concerns the predictability, the safety, the management of illegal activities, as well as the environmental issues that are raised. In this regard, maritime surveillance contributes the most to the discovery of patterns, the decoding of which is a valuable tool for a variety of actions and actions, such as collision avoidance, traffic management and the reduction of marine pollution, among others. Among the various surveillance tools widely used in recent years, the Automatic Identification System (AIS) has established itself as one of the most reliable surveillance technologies found on most ship types. This data consists of kinematic and static messages, providing extremely useful information regarding the movement and identity of ships. In the context of this work, methods and algorithms are proposed to achieve the objectives of processing and analyzing marine data, such as the aforementioned cases. In particular, two spatio-clustering algorithms are presented, which help to locate anchorages and characterize them as long-term or short-term based on their duration (Spatio-Temporal Clustering – STC) and as statistically significant or not (hot/cold spots), taking into account the number of ships for a certain period of time. Finally, experiments were performed using real AIS data received from the antenna of the University of Piraeus in the year 2018 and spatially concern the wider area of the Saronic Gulf. In addition, simulated data were used to conduct escalation experiments.

## Overview
In the scope of the dissertation, two methods were implimented that detect spatial clusters (anchorages) and characterize them either as **long/short-term** based on their temporal persistence `Spatio-Temporal Clustering (STC)`, or as **hot/cold spot** based on their spatio-temporal autocorrelation `Spatio-Temporal Hot Spots (STHS)`. Both methods got inspired by the grid-based method `Trajectory Hot Spots`[^1], where cells in a spatio-temporal partition are clustered and statistically characterized based on the total duration of vessel passing. Anchorages invovles static 

### Preprocessing

### Spatio-Temporal Clustering (STC)

### Spatio-Temporal Hot Spots (STHS)

[^1]: P. Nikitopoulos, A. Paraskevopoulos, C. Doulkeridis, N. Pelekis and Y. Theodoridis, "Hot Spot Analysis over Big Trajectory Data," 2018 IEEE International Conference on Big Data (Big Data), 2018, pp. 761-770, doi: 10.1109/BigData.2018.8622376. 


## Keywords
Maritime Surveillance, Clustering, Statistically Significant Regions, Anchorages, AIS, Spatio-temporal Dimensions

## Quick usage
- Create a conda enviroment using the `environment.yml`
- Execute all cells in notebooks `test_ais.ipynb` and `test_dummy.ipynb`

## README under construction...
