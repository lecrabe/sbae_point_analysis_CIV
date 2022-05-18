import ee
import json
import time
import pandas as pd
import geopandas as gpd
from pathlib import Path
from godale import Executor
from datetime import timedelta

from helpers.ee.get_time_series import get_time_series
from helpers.ee.util import processing_grid
from helpers.ee.landsat.landsat_collection import landsat_collection
from helpers.ee.ccdc import extract_ccdc
from helpers.ee.global_products import sample_global_products_cell

from helpers.ts_analysis.cusum import run_cusum_deforest, cusum_deforest
from helpers.ts_analysis.bfast_wrapper import run_bfast_monitor
from helpers.ts_analysis.bootstrap_slope import run_bs_slope
from helpers.ts_analysis.timescan import run_timescan_metrics
from helpers.ts_analysis.helpers import subset_ts

def get_change_data(aoi, fc, config_dict):
    
    outdir = Path(config_dict['work_dir'])
    if outdir is None:
        outdir = Path.home().joinpath('module_results/sbae_point_analysis')
    
    outdir.mkdir(parents=True, exist_ok=True)
    config_file = str(outdir.joinpath("config.json"))
    
    with open(config_file, "w") as f:
        json.dump(config_dict, f)

    config_file = outdir.joinpath('config.json')
    grid_file = Path(config_file.parent).joinpath('.grid')
    
    if grid_file.exists():
        with open(grid_file, 'r') as f:
            old_grid = f.read()
        
        if config_dict['ts_params']['grid_size'] != float(old_grid):
            print(' You changed grid size. For this all temporary results are now being deleted.')
            for file in Path(config_dict['work_dir']).glob('*tmp'):
                file.unlink()
    else:
        with open(grid_file, 'w') as f:
            f.write(str(config_dict['ts_params']['grid_size']))
    
    # create image collection (not being changed)
    band = config_dict['ts_params']['band']
    lsat = landsat_collection(
        config_dict['ts_params']['start_date'], 
        config_dict['ts_params']['end_date'], 
        aoi, 
        bands=['green', 'red', 'nir', 'swir1', 'swir2', band]
    )
    
    # get parameters from configuration file
    ts_params = config_dict['ts_params']
    band = ts_params['band']
    sat = ts_params['satellite']
    start_hist = ts_params['start_date']
    start_mon = ts_params['start_monitor']
    end_mon = ts_params['end_date']
    grid_size = ts_params['grid_size']

    # get algorithms from config file
    bfast = config_dict['bfast_params']['run']
    cusum = config_dict['cusum_params']['run']
    ccdc = config_dict['ccdc_params']['run']
    ts_metrics = config_dict['ts_metrics_params']['run']
    bs_slope = config_dict['bs_slope_params']['run']
    glb_prd = config_dict['global_products']['run']

    # create namespace for tmp and outfiles
    param_string = f'{sat}_{band}_{start_hist}_{start_mon}_{end_mon}_{grid_size}'
        
    def cell_computation(args):
        
        idx, cell, config_file = args
        with open(config_file, "r") as f:
            config_dict = json.load(f)
        
        
        # get parameters from configuration file
        ts_params = config_dict['ts_params']
        band = ts_params['band']
        sat = ts_params['satellite']
        start_hist = ts_params['start_date']
        start_mon = ts_params['start_monitor']
        end_mon = ts_params['end_date']
        grid_size = ts_params['grid_size']
        
        # get algorithms from config file
        bfast = config_dict['bfast_params']['run']
        cusum = config_dict['cusum_params']['run']
        ccdc = config_dict['ccdc_params']['run']
        ts_metrics = config_dict['ts_metrics_params']['run']
        bs_slope = config_dict['bs_slope_params']['run']
        glb_prd = config_dict['global_products']['run']
        
        # create namespace for tmp and outfiles
        param_string = f'{sat}_{band}_{start_hist}_{start_mon}_{end_mon}_{grid_size}'
        tmp_file = outdir.joinpath(f'tmp_results_{idx}_{param_string}.pickle')
        tmp_empty_file = outdir.joinpath(f'tmp_noresults_{idx}_{param_string}.txt')
        
        # check if already been calculated
        if tmp_file.exists() or tmp_empty_file.exists():
            print(f' Grid cell {idx} already has been extracted. Going on with next grid cell.')    
            return
        
        # get start time
        start_time = time.time()
        
        df = None
        point_id_name = config_dict['ts_params']['point_id']
        
        # get geometry of grid cell and filter points for that
        nr_of_points = fc.filterBounds(cell).size().getInfo()
        if nr_of_points > 0:
            print(f' Processing gridcell {idx}')            

            # get the timeseries data
            if bfast or cusum or ts_metrics or bs_slope:
            
                # extract time-series
                df = get_time_series(lsat.select(band), fc, cell, config_dict)
                
                # run bfast
                df = run_bfast_monitor(df, config_dict['bfast_params']) if bfast else df

                ### THINGS WE RUN WITHOUT HISTORIC PERIOD #####

                # we cut ts data to monitoring period only
                df[['dates', 'ts', 'mon_images']] = df.apply(
                    lambda row: subset_ts(row, config_dict['ts_params']['start_monitor']), axis=1, result_type='expand'
                )
                
                # run cusum
                df = run_cusum_deforest(df, config_dict['cusum_params']) if cusum else df
                
                # run timescan metrics
                df = run_timescan_metrics(df, config_dict['ts_metrics_params']) if ts_metrics else df

                # run bs_slope
                df = run_bs_slope(df, config_dict['bs_slope_params']) if bs_slope else df
            
            # run ccdc
            if ccdc:
                ccdc_df = extract_ccdc(lsat, fc, cell, config_dict)
                if df is not None:
                    df = pd.merge(
                        ccdc_df[['point_id', 'ccdc_change_date', 'ccdc_magnitude']], 
                        df, 
                        on='point_id'
                    )
                else:
                    df = ccdc_df[['point_id', 'ccdc_change_date', 'ccdc_magnitude', 'geometry']]

            if glb_prd:
                glo_products_df = sample_global_products_cell(fc, cell, config_dict)
                
                if df is not None:
                    df = pd.merge(glo_products_df[[
                         'point_id', 'esa_lc20','gfc_tc00','gfc_loss','gfc_year','gfc_gain',
                         'tmf_main_cl','tmf_subtypes','tmf_1990','tmf_1995','tmf_2000',
                         'tmf_2005','tmf_2010','tmf_2015','tmf_2020','tmf_def_yr','tmf_deg_yr']],
                          df, 
                          on='point_id'
                    )
                else:
                    df = glo_products_df[['LON', 'LAT', 'PLOTID',
                         'esa_lc20','gfc_tc00','gfc_loss','gfc_year','gfc_gain',
                         'tmf_main_cl','tmf_subtypes','tmf_1990','tmf_1995','tmf_2000',
                         'tmf_2005','tmf_2010','tmf_2015','tmf_2020','tmf_def_yr','tmf_deg_yr',
                         'geometry']]
            
            #write to tmp pickle file
            df.to_pickle(tmp_file)

            # stop timer and print runtime
            elapsed = time.time() - start_time
            print(f' Grid cell {idx} with {nr_of_points} points done in: {timedelta(seconds=elapsed)}')    
            
        elif nr_of_points == 0:
            with open(tmp_empty_file, 'w') as f:
                f.write('0 points')
            print(f' Grid cell {idx} does not contain any points. Going on with next grid cell.')    
        elif nr_of_points == -1:
            with open(tmp_empty_file, 'w') as f:
                f.write('0 points')  
            print(f' No point data could been extracted from grid cell {idx}. Going on with next grid cell.')        
    
    # create a grid
    grid, grid_fc = processing_grid(aoi, config_dict['ts_params']['grid_size'], config_dict['ts_params']['grid_size'])
    print(f' Parallelizing time-series extraction on {str(config_dict["workers"])} threads for a total of {len(grid)} grid cells.')
    
    args_list = [(*l, config_file) for l in list(enumerate(grid))]
    
    # ---------------debug line--------------------------
    #cell_computation([14, grid[14], config_file])
    # ---------------debug line end--------------------------
    
    executor = Executor(executor="concurrent_threads", max_workers=config_dict["workers"])
    for i, task in enumerate(executor.as_completed(
        func=cell_computation,
        iterable=args_list
    )):
        try:
            task.result()
        except ValueError:
            print("gridcell task failed")
    
    files = list(outdir.glob('tmp_results*.pickle'))
    gdf = pd.read_pickle(files[0])
    df = pd.DataFrame(columns=gdf.columns)
    
    # collect all data into a single dataframe
    for file in outdir.glob('tmp_results*.pickle'):
        df2 = pd.read_pickle(file)
        df = pd.concat([df, df2], ignore_index=True)
    
    # try to turn columsn into numerical, where possible
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except:
            pass
    
    
    out_gpkg = outdir.joinpath(f'results_{param_string}.gpkg')
    out_pckl = outdir.joinpath(f'results_{param_string}.pickle')
    # write to pickle with all ts and dates
    df.to_pickle(out_pckl)
    
    # write to geo file
    if 'dates' in df.columns:
        gdf = gpd.GeoDataFrame(
                    df.drop(['dates', 'ts'], axis=1), 
                    crs="EPSG:4326", 
                    geometry=df['geometry']
                )
    else:
        gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=df['geometry'])
        
    ## write to output and return df# write to output and return df
    gdf.to_file(out_gpkg, driver='GPKG')
    
    
    for file in Path(config_dict['work_dir']).glob('tmp*'):
        file.unlink()
    print(" Processing has been finished successfully. Check for final_results files in your output directory.")