from pathlib import Path 

import ee
import geopandas as gpd

def rename_TMF(band):
    return ee.String(band).replace("Dec", "tmf_",'g')


def sample_global_products(points_fc, outfile):

    ## Global Forest Change (Hansen et al., 2013)
    gfc_col = ee.Image(
        'UMD/hansen/global_forest_change_2020_v1_8'
    ).select(
        ['treecover2000','loss','lossyear','gain'],
        ['gfc_tc00','gfc_loss','gfc_year','gfc_gain']
    )

    ## ESA WorldCover 2020
    esa_20  = ee.Image('ESA/WorldCover/v100/2020').rename('esa_lc20')

    ## Tropical Moist Forest - JRC 2021
    tmf_annual= ee.ImageCollection('projects/JRC/TMF/v1_2020/AnnualChanges').mosaic()
    tmf_ann_n = tmf_annual.rename(tmf_annual.bandNames().map(rename_TMF))

    tmf_subtp = ee.ImageCollection('projects/JRC/TMF/v1_2020/TransitionMap_Subtypes').mosaic().rename('tmf_subtypes')
    tmf_main  = ee.ImageCollection('projects/JRC/TMF/v1_2020/TransitionMap_MainClasses').mosaic().rename('tmf_main_cl')
    tmf_deg   = ee.ImageCollection('projects/JRC/TMF/v1_2020/DegradationYear').mosaic().rename('tmf_deg_yr')
    tmf_def   = ee.ImageCollection('projects/JRC/TMF/v1_2020/DeforestationYear').mosaic().rename('tmf_def_yr')

    ##  COMBINE COLLECTIONS
    glo_ds = (esa_20
            .addBands(gfc_col)
            .addBands(tmf_subtp)
            .addBands(tmf_main)
            .addBands(tmf_deg)
            .addBands(tmf_def)
            .addBands(tmf_ann_n)
        )

    ## EXTRACT THE DATA TO THE POINTS
    sampled_points = glo_ds.reduceRegions(**{
      "reducer": ee.Reducer.first(),
      "collection": points_fc
    }).select(
        ['point_id','esa_lc20','gfc_tc00','gfc_loss','gfc_year','gfc_gain',
         'tmf_main_cl','tmf_subtypes','tmf_1990','tmf_1995','tmf_2000',
         'tmf_2005','tmf_2010','tmf_2015','tmf_2020','tmf_def_yr','tmf_deg_yr',
         '.geo']);
    
    ## MAKE AS A DATAFRAME AND EXPORT
    json = sampled_points.getInfo()
    gdf = gpd.GeoDataFrame.from_features(json)
    gdf['LON'] = gdf['geometry'].x
    gdf['LAT'] = gdf['geometry'].y
    
    # sort columns for CEO output
    gdf['PLOTID'] = gdf['point_id']
    cols = gdf.columns.tolist()
    cols = [e for e in cols if e not in ('LON', 'LAT', 'PLOTID')]
    new_cols = ['LON', 'LAT', 'PLOTID'] + cols
    gdf = gdf[new_cols]
    gdf.to_csv(Path.home()/outfile, index=False)
    return gdf