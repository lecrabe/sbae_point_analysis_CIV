import ee
import numpy as np
import pandas as pd
from retry import retry
from godale import Executor


fake_geometry = ee.Geometry.Polygon([
    [[
        14.642002058678267,
        -1.5396013709395364
      ],
      [
        21.497470808678266,
        -1.5396013709395364
      ],
      [
        21.497470808678266,
        2.326761024980699
      ],
      [
        14.642002058678267,
        2.326761024980699
      ],
      [
        14.642002058678267,
        -1.5396013709395364
      ]
    ]
  ])

def zip2Image(element):
  
    value = ee.List(element).get(0)
    year = ee.List(element).get(1)
    
    return (ee.Image.constant(ee.Number.parse(value))
      .rename('ndvi')
      .clip(fake_geometry)
      .set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis())
      .toFloat()
           )

@retry(tries=5, delay=1, backoff=2)
def extract_landtrendr(args_list):
    
    years_list, ts_list, point_id, landtrendr_params = args_list
    ts = ee.List([str(value) for value in ts_list])
    dates = ee.List([int(year) for year in years_list])

    ts = ts.zip(dates)
    startYear = np.min(years_list)
    endYear = np.max(years_list)

    tsee = ee.ImageCollection(
      ts.map(zip2Image)
    )

    
    landtrendr_params.update(timeSeries=tsee)
    landtrendr_params.pop('run', None)

    lt = ee.Algorithms.TemporalSegmentation.LandTrendr(**landtrendr_params).select(["LandTrendr"])

    vertexMask = lt.arraySlice(0, 3, 4); # slice out the 'Is Vertex' row - yes(1)/no(0)
    vertices = lt.arrayMask(vertexMask); # use the 'Is Vertex' row as a mask for all rows

    left = vertices.arraySlice(1, 0, -1);    # slice out the vertices as the start of segments
    right = vertices.arraySlice(1, 1, None); # slice out the vertices as the end of segments
    startYear = left.arraySlice(0, 0, 1);    # get year dimension of LT data from the segment start vertices
    startVal = left.arraySlice(0, 2, 3);     # get spectral index dimension of LT data from the segment start vertices
    endYear = right.arraySlice(0, 0, 1);     # get year dimension of LT data from the segment end vertices 
    endVal = right.arraySlice(0, 2, 3);      # get spectral index dimension of LT data from the segment end vertices

    dur = endYear.subtract(startYear);       # subtract the segment start year from the segment end year to calculate the duration of segments 
    mag = endVal.subtract(startVal);         # substract the segment start index value from the segment end index value to calculate the delta of segments
    rate = mag.divide(dur);                  # calculate the rate of spectral change

    segInfo = (
        ee.Image.cat([startYear.add(1), endYear, startVal, endVal, mag, dur, rate])
            .toArray(0)
            .mask(vertexMask.mask())
    )

    distDir = -1;

    sortByThis = segInfo.arraySlice(0,4,5).toArray(0).multiply(-1); # need to flip the delta here, since arraySort is working by ascending order
    segInfoSorted = segInfo.arraySort(sortByThis); # sort the array by magnitude
    bigDelta = segInfoSorted.arraySlice(1, 0, 1); # get the first segment in the sorted array (greatest magnitude vegetation loss segment)

    bigDeltaImg = ee.Image.cat(bigDelta.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['yod']]),
    bigDelta.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['endYr']]),
    bigDelta.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['startVal']]).multiply(distDir),
    bigDelta.arraySlice(0,3,4).arrayProject([1]).arrayFlatten([['endVal']]).multiply(distDir),
    bigDelta.arraySlice(0,4,5).arrayProject([1]).arrayFlatten([['mag']]).multiply(distDir),
    bigDelta.arraySlice(0,5,6).arrayProject([1]).arrayFlatten([['dur']]),
    bigDelta.arraySlice(0,6,7).arrayProject([1]).arrayFlatten([['rate']]).multiply(distDir));

    distMask =  bigDeltaImg.select(['mag']).lt(1000).And(bigDeltaImg.select(['dur']).lt(5));

    bigFastDist = bigDeltaImg  #.mask(distMask).int16(); // need to set as int16 bit to use connectedPixelCount for minimum mapping unit filter

    landtrendr = bigFastDist.select(['mag', 'dur', 'yod', 'rate', 'endYr']).clip(fake_geometry).reduceRegion(**{
      'reducer': ee.Reducer.first(),
      'scale': 3000
    }).getInfo()

    return landtrendr['mag'], landtrendr['dur'], landtrendr['yod'], landtrendr['rate'],  landtrendr['endYr'], point_id


def run_landtrendr(df, landtrendr_params):

    args_list, d = [], {}

    for i, row in df.iterrows():
        
        # get years of ts
        years = np.unique([date.year for date in row.dates])
        
        # get mean value for each year
        ts_yearly = []
        for year in years:

            idx = np.array([True if date.year == year else False for date in row.dates])
            ts_yearly.append(np.nanmean(np.array(row.ts)[idx]))

        args_list.append([years, ts_yearly, row.point_id, landtrendr_params])
        
        
    executor = Executor(executor="concurrent_threads", max_workers=16)
    for i, task in enumerate(executor.as_completed(
        func=extract_landtrendr,
        iterable=args_list
    )):
        try:
            d[i] = list(task.result())
        except ValueError:
            print("landtrendr task failed")

    landtrendr_df = pd.DataFrame.from_dict(d, orient='index')
    landtrendr_df.columns = ['ltr_magnitude', 'ltr_dur', 'ltr_yod', 'ltr_rate', 'ltr_end_year', 'point_id']
    return pd.merge(df, landtrendr_df, on='point_id')
