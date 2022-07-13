import ee
import requests
import numpy as np
import pandas as pd
import geopandas as gpd
from retry import retry


def extract_landtrendr(time_series, config_dict):
    
    ts = ee.List([str(value) for value in ts_list])
    dates = ee.List([int(year) for year in years_list])

    ts = ts.zip(dates)
    startYear = np.min(years_list)
    endYear = np.max(years_list)

    tsee = ee.ImageCollection(
      ts.map(zip2Image)
    )

    runParams = { 
        'maxSegments':            6,
        'spikeThreshold':         0.9,
        'vertexCountOvershoot':   3,
        'preventOneYearRecovery': True,
        'recoveryThreshold':      0.25,
        'pvalThreshold':          0.05,
        'bestModelProportion':    0.75,
        'minObservationsNeeded':  3,
        'timeSeries':            tsee
    }

    lt = ee.Algorithms.TemporalSegmentation.LandTrendr(**runParams).select(["LandTrendr"])

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

    return bigFastDist.select('mag').clip(geometry).reduceRegion(**{
      'reducer': ee.Reducer.first(),
      'scale': 3000
    }).getInfo()