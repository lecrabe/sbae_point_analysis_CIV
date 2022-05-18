from pathlib import Path
import ee
import geemap
import time
import pandas as pd
import geopandas as gpd
from shapely.geometry import box, Point
import numpy as np


def random_point(geometry):
    bounds = geometry.bounds
    x = (bounds[2] - bounds[0]) * np.random.random_sample(1) + bounds[0]
    y = (bounds[3] - bounds[1]) * np.random.random_sample(1) + bounds[1]
    return Point(x, y)
    
    
def generate_samples(aoi, spacing, crs='ESRI:54017', sampling_strategy='systematic'):

    if isinstance(aoi, ee.FeatureCollection):
        aoi = geemap.ee_to_geopandas(aoi)
    
    # reproject
    aoi = aoi.dissolve().to_crs(crs)
    aoi_geom = aoi.iloc[0]['geometry']
    
    # get bounds
    bounds = aoi.bounds
    
    # get orgiin point
    originx = bounds.minx.values[0]
    originy = bounds.miny.values[0]

    # get widht and height of aoi bounds
    width = bounds.maxx - bounds.minx 
    height = bounds.maxy - bounds.miny

    # calculate how many cols and row are those
    columns = int(np.floor(float(width) / spacing))
    rows = int(np.floor(float(height) / spacing))
    
    # create grid cells
    print("Creating grid cells")
    i, l = 1, []
    for column in range(0, columns + 1):
        x = originx + (column * spacing)
        for row in range(0, rows + 1):
            y = originy + (row * spacing)
            cell = box(x, y, x+spacing, y+spacing)
            if cell.intersects(aoi_geom):
                l.append(cell)
                i += 1
    
    # and turn into geodataframe
    print("Turning grid cells into GeoDataFrame...")
    df = pd.DataFrame(l, columns=['geometry'])
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=crs) 
    
    # add points
    print("Creating sampling points...")
    if sampling_strategy == 'systematic':
        # take centroid
        gdf['sample_points'] = gdf.geometry.centroid
        
    elif sampling_strategy == 'random': 
        # create rand points in each grid
        gdf['sample_points'] = gdf.geometry.apply(lambda shp: random_point(shp))
        
    # add point id
    print("Adding a unique point ID...")
    gdf['point_id'] = [i for i in range(len(gdf.index))]
    
    # divide to grid and point df
    grid_gdf = gdf.drop(['sample_points'], axis=1)
    gdf['geometry'] = gdf['sample_points']
    point_gdf = gdf.drop(['sample_points'], axis=1)
    
    print('Remove points outside AOI...')
    point_gdf = point_gdf[point_gdf.geometry.within(aoi_geom)]
    return grid_gdf, point_gdf


def split_dataframe(df, chunk_size = 25000): 
        chunks = list()
        num_chunks = len(df) // chunk_size + 1
        for i in range(num_chunks):
            chunks.append(df[i*chunk_size:(i+1)*chunk_size])
        return chunks
    

def upload_to_ee(gdf, asset_name):
    
    # get users asset root
    asset_root = ee.data.getAssetRoots()[0]['id']

    if len(gdf) > 25000:
        print('Need to run splitted upload routine as dataframe has more than 25000 rows')

        try:
            # create temporary folder
            ee.data.createAsset({'type': 'folder'}, f'{asset_root}/tmp_sbae')
        except:
            pass

        # upload chunks of data to avoi max upload
        chunks = split_dataframe(point_df)
        tasks = []

        for i, chunk in enumerate(chunks):

            point_fc = geemap.geopandas_to_ee(chunk.to_crs("EPSG:4326"))
            exportTask = ee.batch.Export.table.toAsset(
                collection = point_fc,
                description = f'sbae_part_{i}',
                assetId = f'{asset_root}/tmp_sbae/points_{i}'
            )

            exportTask.start()
            tasks.append(exportTask)

        # cha on status
        
        finished=False
        while finished == False:
            time.sleep(30)
            for task in tasks:
                state = task.status()['state']
                if state == 'COMPLETED':
                    finished = True
                else:
                    finished = False
                    break

        # merge assets
        print('aggregate to final')
        child_assets = ee.data.listAssets({'parent': f'{asset_root}/tmp_sbae'})['assets']
        for i, ass in enumerate(child_assets):
            if i == 0:
                point_fc = ee.FeatureCollection(ass['id'])
            else:
                point_fc = point_fc.merge(ee.FeatureCollection(ass['id']))

        print('export final')
        # export to final
        exportTask = ee.batch.Export.table.toAsset(
                collection = point_fc,
                description = f'sbae_aggregated_table',
                assetId = f'{asset_root}/{asset_name}'
            )

        exportTask.start()

        finished=False
        while finished == False:
            time.sleep(30)
            state = exportTask.status()['state']

            if state == 'COMPLETED':
                finished = True
            else:
                finished = False

        print('delete temporary assets')
        child_assets = ee.data.listAssets({'parent': f'{asset_root}/tmp_sbae'})['assets']
        for i, ass in enumerate(child_assets):
            ee.data.deleteAsset(ass['id'])  

        ee.data.deleteAsset(f'{asset_root}/tmp_sbae')
        print(f' Upload completed. You can find the samples at {asset_name}')
        
    else:
        
        # turn into FC
        point_fc = geemap.geopandas_to_ee(gdf.to_crs("EPSG:4326"))
        
        print(' Exporting to asset')
        # export to final
        exportTask = ee.batch.Export.table.toAsset(
                collection = point_fc,
                description = f'sbae_aggregated_table',
                assetId = f'{asset_root}/{asset_name}'
            )

        exportTask.start()

        finished=False
        while finished == False:
            time.sleep(30)
            state = exportTask.status()['state']

            if state == 'COMPLETED':
                finished = True
            else:
                finished = False

                
def save_locally(gdf, ceo_csv=True, gpkg=True, outdir=None):
    
    if not outdir:
        outdir = Path.home().joinpath('module_results/sbae_point_analysis')
        
    outdir.mkdir(parents=True, exist_ok=True)
    
    print(f' Saving outputs to {outdir}')
    gdf['LON'] = gdf['geometry'].x
    gdf['LAT'] = gdf['geometry'].y
    
    # sort columns for CEO output
    gdf['PLOTID'] = gdf['point_id']
    cols = gdf.columns.tolist()
    cols = [e for e in cols if e not in ('LON', 'LAT', 'PLOTID')]
    new_cols = ['LON', 'LAT', 'PLOTID'] + cols
    gdf = gdf[new_cols]
    
    if ceo_csv:
        gdf.to_csv(outdir.joinpath('01_sbae_points.csv'), index=False)
        
    if gpkg:
        gdf.to_file(outdir.joinpath('01_sbae_points.gpkg'), driver='GPKG')
        