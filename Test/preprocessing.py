import os
import geopandas as gpd
from sqlalchemy.engine import create_engine
from shapely.geometry import LineString, Point
import json

def _get_engine():
    '''Create a Postgres engine based on dbcon.json DB credentials.
    '''
    data = json.load(open(os.path.join(os.getcwd(), "dbcon.json"), "r"))
    username = data['username']
    pswd = data['pswd']
    host = data['host']
    port = data['port']
    db = data['db']
    return create_engine(f"postgresql://{username}:{pswd}@{host}:{port}/{db}")


def _read_db():
    '''Reads the DB using the query in dbcon.json.
    '''
    query = json.load(open(os.path.join(os.getcwd(), "dbcon.json"), "r"))['query']
    gdf = gpd.read_postgis(query, _get_engine(), geom_col='geom', crs=4326)
    gdf = gdf.sort_values(['mmsi', 't']).reset_index(drop=True)
    return gdf


def _apply_find_gaps(df, dist_threshold):
    '''Returns the gap indexes of a trajectory, based on dist_threshold
    '''
    gaps = df.loc[df['geom'].distance(df['geom'].shift()).values / 1000 > dist_threshold].index.to_list()

    if len(df) > 1 and df.iloc[1]['geom'].distance(df.iloc[0]['geom']) / 1000 > dist_threshold:
        gaps.extend([df.index.min()])
    gaps.extend([df.index.min(), df.index.max()])
    gaps.sort()
    return gaps


def _apply_split_to_trajectories(object_df, gaps):
    '''Splits a trajectory based on the gaps.
    '''
    if len(object_df) == 1:
        return None
    if len(object_df)<3:
        object_df.loc[:, 'tid'] = 0
        return object_df

    for i, (start_ind, end_ind) in enumerate(list(zip(gaps, gaps[1:]))):
        object_df.loc[start_ind:end_ind, 'tid'] = i

    return object_df


def _apply_create_lines(gdf):
    '''Creates a Listestring object, given a trajectory in GeoDataFrame.
    '''
    
    start_t = min(gdf['t'])
    end_t = max(gdf['t'])
    if len(gdf.geom.tolist()) <= 1:
        return None
    else:
        line = LineString(gdf.geom.tolist())
        return gpd.GeoDataFrame({"t_s": [start_t], "t_e": [end_t], "geom": [line]}, geometry="geom", crs=4326)   


def _apply_calc_m(line):
    '''Calculates the midpoint of a Line.
    '''

    x_m = sum(x for x in line.geom.iloc[0].xy[0]) / len(line.geom.iloc[0].coords)
    y_m = sum(y for y in line.geom.iloc[0].xy[1]) / len(line.geom.iloc[0].coords)
    m = Point(x_m, y_m)
    return gpd.GeoDataFrame({"t_s": [line['t_s'].iloc[0]], "t_e": [line['t_e'].iloc[0]], "geom": [m]}, geometry="geom", crs=4326)


def preprocessing():
    '''Loads the data from the DB and returns the cleaned data, the created trajectories and the midpoints
    '''
    gdf = _read_db()

    dist_threshold = 10 #km
    gdf.to_crs(2100, inplace=True)
    gaps = gdf.groupby('mmsi').apply(lambda grp: _apply_find_gaps(grp, dist_threshold))
    gdf_clean = gdf.groupby('mmsi').apply(lambda grp: _apply_split_to_trajectories(grp, gaps.loc[int(grp.name)]))
    gdf_clean.set_geometry('geom', inplace=True)
    gdf_clean.to_crs(4326, inplace=True)
    gdf_clean.reset_index(drop=True,inplace=True)

    lines_df = gdf_clean.groupby(['mmsi', 'tid']).apply(lambda g: _apply_create_lines(g)).dropna().reset_index().drop("level_2", axis=1)
    
    m_df = lines_df.groupby(['mmsi', 'tid']).apply(lambda line: _apply_calc_m(line)).reset_index().drop('level_2', axis=1)

    m_df.to_csv('./ms.csv', index=False)

    return gdf_clean, lines_df, m_df