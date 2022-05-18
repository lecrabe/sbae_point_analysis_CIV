import numpy as np
from datetime import datetime as dt
import pandas as pd
import seaborn as sns


def subset_ts(row, start_monitor):
    """ Helper function to extract only monitoring period
    """
    
    # create index for monitoring period
    idx = row.dates > dt.strptime(start_monitor, '%Y-%m-%d')
    
    # subset dates
    dates = row.dates[idx]
    
    # subset ts data
    ts = np.array(row.ts)[idx].tolist()
    
    # get new image length
    images = len(ts)
    
    return dates, ts, images


def plot_timeseries(pickle_file, point_id, point_id_name='point_id'):
    
    df = pd.read_pickle(pickle_file)
    dates = df[df[point_id_name] == point_id].dates.values[0]
    ts = np.array(df[df[point_id_name] == point_id].ts.values[0])

    sns.scatterplot(x=dates, y=ts)