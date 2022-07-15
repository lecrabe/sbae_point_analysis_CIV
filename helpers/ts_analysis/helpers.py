import numpy as np
from scipy import stats
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


def rolling_mean(dates, ts, interval='60d'):
    tmp_df = pd.DataFrame(data=ts, index=pd.DatetimeIndex(dates), columns=['ts'])
    return tmp_df.rolling(interval).mean().ts.tolist()


def smooth_ts(df):

    df['ts'] = df.apply(lambda x: rolling_mean(x.dates, x.ts), axis=1)
    return df


def outlier_removal(dates, ts):
    ts = np.array(ts).astype(float)
    z_score = np.abs(stats.zscore(ts, axis=0))
    ts[z_score > 3] = np.nan
    tmp_df = pd.DataFrame(data=ts, index=pd.DatetimeIndex(dates), columns=['ts'])
    tmp_df = tmp_df.dropna()
    return tmp_df.index, tmp_df.ts.tolist()


def remove_outliers(df):

    df[['dates', 'ts']] = df.apply(lambda x: outlier_removal(x.dates, x.ts), axis=1, result_type='expand')
    return df


def plot_timeseries(pickle_file, point_id, point_id_name='point_id'):
    
    df = pd.read_pickle(pickle_file)
    dates = df[df[point_id_name] == point_id].dates.values[0]
    ts = np.array(df[df[point_id_name] == point_id].ts.values[0])

    sns.scatterplot(x=dates, y=ts)