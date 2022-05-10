"""
This is modified python script for calculating (and saving) visualizations from `vis.py` and `vis_dask.py`
"""
import rvt.vis
import rvt.vis_dask
import rvt.blend
import rvt.default
import os
import dask
from dask.distributed import Client, LocalCluster, Lock
import time
import json
import psutil
import logging
logger = logging.getLogger(__name__)
#import dask_memusage
import gc

input_dir_path = "input_dir/large/"
output_dir_path = "output_dir/large/"
# chunksize_mib = ["2MiB", "4MiB", "6MiB", "8MiB", "12MiB", "16MiB", "32MiB" ]
n_dir_list = [1, 4, 8, 16]
chunksize_mib = ["32MiB","64MiB", "128MiB" , "256MiB", "512MiB"]

def get_numpy_result(input_dir_path,n_dir):
    dem_list = os.listdir(input_dir_path)
    for input_dem_name in dem_list:
        if ".tif" not in input_dem_name:  # preskoči če se file ne konča .tif
            continue
        input_dem_path = os.path.join(input_dir_path, input_dem_name)
        out_name = "{}_svf_ndir{}.tif".format(input_dem_name.rstrip(".tif"), n_dir)
        save_render_path= os.path.abspath(os.path.join(output_dir_path, out_name))
    
        input_dict_np_arr = rvt.default.get_raster_arr(input_dem_path)
        input_np_arr = input_dict_np_arr ['array']
        x_res = input_dict_np_arr['resolution'][0]      
        y_res = input_dict_np_arr['resolution'][1] 
        no_data = input_dict_np_arr['no_data']
        ve_factor = 1
        dict_svf_asvf_opns = rvt.vis.sky_view_factor(dem=input_np_arr, resolution=x_res, compute_svf=True,
                                                compute_asvf=False, compute_opns=False,
                                                svf_n_dir=n_dir, svf_r_max=10,
                                                svf_noise=0, asvf_dir=315,
                                                asvf_level=1, ve_factor=ve_factor,
                                                no_data = no_data)
        arr, = dict_svf_asvf_opns.values() 
        # rvt.default.save_raster(src_raster_path= input_dem_path, out_raster_path=save_render_path,
        #                                 out_raster_arr=arr, no_data=None, e_type=6)
    return arr
        

def get_dask_result(input_dir_path, test_chunksize, workr, thrds, n_dir): 
    dem_list = os.listdir(input_dir_path)
    for input_dem_name in dem_list:
        if ".tif" not in input_dem_name:  # preskoči če se file ne konča .tif
            continue
        input_dem_path = os.path.join(input_dir_path, input_dem_name)
        out_name = "{}_svf_ndir{}_chunk{}_w{}_t{}.tif".format(input_dem_name.rstrip(".tif"), n_dir, test_chunksize, workr, thrds)
        save_render_path= os.path.abspath(os.path.join(output_dir_path, out_name))
    
        input_dict_da_arr = rvt.default.get_raster_dask_arr(input_dem_path)
        input_da_arr = input_dict_da_arr['array']                                   ## compute entire array
        # h, w = input_dict_da_arr['array'].shape
        # da_arr = input_dict_da_arr['array'][0 : h//4, w//4: 3 * w//4]           ## or compute only desired region of the array, e.g. part of very large raster

        x_res = input_dict_da_arr['resolution'][0]      
        y_res = input_dict_da_arr['resolution'][1] 
        no_data = input_dict_da_arr['no_data']
        ve_factor = 1
        arr = rvt.vis_dask.dask_sky_view_factor(input_dem=input_da_arr, resolution=x_res,compute_svf=True,
                                                    compute_asvf=False, compute_opns=False,
                                                    svf_n_dir= n_dir, svf_r_max=10,
                                                    svf_noise=0, asvf_dir=315,
                                                    asvf_level=1, ve_factor=ve_factor,
                                                    no_data = no_data)
        ##eg svf                                        
        # arr = rvt.vis_dask.dask_slope_aspect(input_dem=input_da_arr, resolution_x=x_res,
        #                                 resolution_y=y_res, ve_factor=ve_factor,
        #                                 output_units='degree')[0,:,:].squeeze()
        rvt.default.dask_save_raster_tif(src_raster_path= input_dem_path, out_raster_path=save_render_path,
                                        out_raster_arr=arr, dtype_to_save= "float32")
    return 0

def txt_output(output_file_path, output_file):
    with open(output_file_path, "w") as outf:
        json.dump(output_file, outf)
    return 0

###============================
# Program start
if __name__ == '__main__':

    ### NUMPY
    # for n in n_dir_list:
    #     print(f' n_directions: {n}')
    #     start = time.time()
    #     try:
    #         get_numpy_result(input_dir_path = input_dir_path, n_dir = n)
    #         end = time.time()
    #         record = {'name': "numpy", 
    #                 'duration': end - start, 
    #                 'n_dir': n}
    #         outfile = os.path.join(output_dir_path, f'n_dir{n}.txt')  ###intermediate results, in case of workers crashing
    #         txt_output(output_file_path=outfile, output_file = record)
    #     except Exception as e:
    #         end = time.time()
    #         logger.exception(f"Exception {e}.")
    #         record = {'name': "numpy", 
    #                 'exception': end - start,  
    #                 'n_dir': n}
    #         outfile = os.path.join(output_dir_path, f'n_dir{n}.txt')  ###intermediate results, in case of workers crashing
    #         txt_output(output_file_path=outfile, output_file = record)
    #         continue
    #     finally:
    #         time.sleep(2)
    #         gc.collect()

    ### DASK
    # cluster = LocalCluster()
    # client = Client(cluster)
    # total_cores = sum(client.ncores().values())            ## = 16
    # client.shutdown()

    tot_memory = psutil.virtual_memory().total >> 30        ## bytes, shift GB (aprox, rounding error)
    total_cores = 16
    workers_threads_pairs = [ (1, 8), (1,16), (2, 4), (2, 8), (4, 12), (4, 4)]
    # workers_threads_pairs = [(1, total_cores), (2, total_cores//2), (4, total_cores//4)]
    # workers_threads_pairs = [(1, total_cores - 2), (2, total_cores//2 - 2), (4, total_cores//4 - 2)] ### don't use all cores, for example - 2 threads per worker

    # dask_memusage.install(cluster.scheduler, f"memusage_{n_workrs}_{chunksize}.csv")
    out=[]
    for n_workrs, n_thrds in workers_threads_pairs:      
        for chunksize in chunksize_mib:
            for n in n_dir_list:
                print(f'workers: {n_workrs}, threads:{n_thrds}, n_directions: {n}, chunksize: {chunksize}')
                # try:
                cluster = LocalCluster(n_workers=n_workrs, threads_per_worker = n_thrds,
                                    # memory_limit = f'{tot_memory//n_workrs}GB')
                                    memory_limit = None)
                # dask_memusage.install(cluster.scheduler, f"memusage_{n_workrs}_{chunksize}.csv")
                client = Client(cluster)
                # client.dashboard_link       ## <- open in browser to see dask diagnostics dashboard live

                start = time.time()
                try:
                    dask.config.set({"array.chunk-size": chunksize})
                    # print(dask.config.get("array.chunk-size"))
                    get_dask_result(input_dir_path = input_dir_path, test_chunksize = chunksize, workr = n_workrs, thrds = n_thrds, n_dir = n)
                    end = time.time()

                    record = {'name': f"n_workers: {n_workrs}, n_threads: {n_thrds}", 
                            'duration': end - start, 
                            'chunk_size': chunksize, 
                            'n_dir': n,
                            'n_workers': n_workrs, 
                            'n_threads': n_thrds}
                    
                    outfile = os.path.join(output_dir_path, f'n_dir{n}_w{n_workrs}_t{n_thrds}_{chunksize}.txt')  ###intermediate results, in case of workers crashing
                    txt_output(output_file_path=outfile, output_file = record)
                
                except Exception as e:
                    end = time.time()
                    logger.exception(f"Exception {e}.")
                    record = {'name': f"n_workers: {n_workrs}, n_threads: {n_thrds}", 
                            'exception': end - start, 
                            'chunk_size': chunksize, 
                            'n_dir': n,
                            'n_workers': n_workrs, 
                            'n_threads': n_thrds}
            
                    outfile = os.path.join(output_dir_path, f'n_dir{n}_w{n_workrs}_t{n_thrds}_{chunksize}.txt')  ###intermediate results, in case of workers crashing
                    txt_output(output_file_path=outfile, output_file = record)
                    continue
                finally:
                    out.append(record)
                    client.shutdown()
                    time.sleep(2)
                    gc.collect()

    outfile = os.path.join(output_dir_path, 'all_results.txt')
    txt_output(output_file_path=outfile, output_file = out)

