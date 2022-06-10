"""
This is modified python script for calculating VAT Combined blender combination with dask distributed, Local Cluster.

Apply ONE blender comb. at a time: general, flat or combined (duration is per chunksize). Flat and general OK individually, combined acting strange for some chunksizes (fix something in create_layer when `image` is not None, is just the order of vat_arr_1 and vat_arr_2 getting mixed up for some reason in this script? If reading VAT_flat and VAT_general from `image_path`, results are ok. Maybe keep just image_path with dask arrays - good practice to save intermediate results to break up workload anyway...)
If saving MULTIPLE combinations: general, flat and combined at the same time (ONLY on dask distributed - cluster) sometimes problem with aquiring locks when writing (random, hard to reproduce). Look into this, preferably without changing the code in `distributed/lock.py` or 
`dask/array/core.py`

"Luminosity" failing tests in case of multiband data! (present at /blender_VAT.json". Must fix.)

(1st layer: VAT general (with opacity set in parameters, blend mode normal), 2nd layer: VAT flat)
Parameters are:
    input_dir_path (dir that contains GeoTIFFs of DEMs)
    output_dir_path (dir where to store calculated GeoTIFFS of VAT combined visualizations)
    general_opacity (opacity of VAT general)
    test_chunksize (chunksize in MiB)
    workr (number of workers)
    thrds (number of threads per worker)
    vat_combination_json_path (path to VAT combination JSON file)
    terrains_sett_json_path (path to terrains settings JSON file)
    nr_processes (number of parallel processes which are calculating and saving output VAT combined)
    save_float (if saves float)
    save_8bit (if saves 8bit)
    save_VAT_general (if saves VAT general)
    save_VAT_flat (if saves VAT flat)
"""
import rvt.vis
import rvt.blend
import rvt.default
import os
import dask
from dask.distributed import Client, LocalCluster, Lock
import time
import json
import psutil
# import dask_memusage
import gc
import logging
logger = logging.getLogger(__name__)

input_dir_path = "input_dir/large/"
output_dir_path = "output_dir/large/"
general_opacity = 50
vat_combination_json_path = "settings/blender_VAT.json"
terrains_sett_json_path = "settings/default_terrains_settings.json"
save_float = True
save_8bit = False
save_VAT_general = True
save_VAT_flat = False
# chunksize_mib = ["8MiB", "16MiB", "32MiB", "64MiB", "128MiB", "256MiB", "512MiB"]  ##try different chunksizes
chunksize_mib = ["16MiB", "32MiB", "64MiB", "128MiB" , "256MiB", "512MiB"] ## larger images
# chunksize_mib = ["1MiB", "2MiB", "4MiB", "8MiB" , "16MiB", "32MiB"]   ##smaller images


def combined_VAT(input_dir_path, output_dir_path, 
                test_chunksize, workr, thrds,          
                general_opacity, vat_combination_json_path=None,
                terrains_sett_json_path=None, save_float=True, save_8bit=False,
                save_VAT_general=False, save_VAT_flat=False):
    if not save_float and not save_8bit:
        raise Exception("save_float and save_8bit are both False!")

    if vat_combination_json_path is None:  # Če ni podan path do vat_comb se smatra da je v settings
        vat_combination_json_path = os.path.abspath(os.path.join("settings", "blender_VAT.json"))

    if terrains_sett_json_path is None:  # Če ni podan path do terrain sett se smatra da je v settings
        terrains_sett_json_path = os.path.abspath(os.path.join("settings",
                                                               "default_terrains_settings.json"))

    # 2 default classes one for VAT general one for VAT flat
    default_1 = rvt.default.DefaultValues()  # VAT general
    default_2 = rvt.default.DefaultValues()  # VAT flat

    # fill no_data and original no data
    default_1.fill_no_data = 0
    default_2.fill_no_data = 0
    default_1.keep_original_no_data = 0
    default_2.keep_original_no_data = 0

    # 2 blender combination classes one for VAT general one for VAT flat
    vat_combination_1 = rvt.blend.BlenderCombination()  # VAT general
    vat_combination_2 = rvt.blend.BlenderCombination()  # VAT flat

    # read combination from JSON
    vat_combination_1.read_from_file(vat_combination_json_path)
    vat_combination_2.read_from_file(vat_combination_json_path)

    # create terrains settings and read it from JSON
    terrains_settings = rvt.blend.TerrainsSettings()
    terrains_settings.read_from_file(terrains_sett_json_path)

    # select single terrain settings for general and flat
    terrain_1 = terrains_settings.select_terrain_settings_by_name("general")  # VAT general
    terrain_2 = terrains_settings.select_terrain_settings_by_name("flat")  # VAT flat

    # apply terrain settings on default and combination
    terrain_1.apply_terrain(default=default_1, combination=vat_combination_1)  # VAT general
    terrain_2.apply_terrain(default=default_2, combination=vat_combination_2)  # VAT flat

    dem_list = os.listdir(input_dir_path)
    for input_dem_name in dem_list:
        if ".tif" not in input_dem_name:  # preskoči če se file ne konča .tif
            continue
        input_dem_path = os.path.join(input_dir_path, input_dem_name)
        out_name = "{}_Archaeological_(VAT_combined)_chunk{}_w{}_t{}.tif".format(input_dem_name.rstrip(".tif"), test_chunksize, workr, thrds)
        out_comb_vat_path = os.path.abspath(os.path.join(output_dir_path, out_name))
        out_comb_vat_8bit_path = out_comb_vat_path.rstrip(".tif") + "_8bit.tif"
        if save_8bit and os.path.isfile(out_comb_vat_8bit_path) and save_float and os.path.isfile(out_comb_vat_path):
            print("{} already exists!".format(out_comb_vat_path))
            print("{} already exists!".format(out_comb_vat_8bit_path))
            continue
        elif save_float and os.path.isfile(out_comb_vat_path) and not save_8bit:  # output VAT comb already exists
            print("{} already exists!".format(out_comb_vat_path))
            continue
        elif save_8bit and os.path.isfile(out_comb_vat_8bit_path) and not save_float:
            print("{} already exists!".format(out_comb_vat_8bit_path))
            continue
        general_combination = vat_combination_1
        flat_combination = vat_combination_2
        general_default = default_1
        flat_default = default_2
        out_comb_vat_general_name = "{}_Archaeological_(VAT_general)_chunk{}_w{}_t{}.tif".format(input_dem_name.rstrip(".tif"), test_chunksize, workr, thrds)
        out_comb_vat_general_path = os.path.abspath(os.path.join(output_dir_path, out_comb_vat_general_name))
        out_comb_vat_flat_name = "{}_Archaeological_(VAT_flat)_chunk{}_w{}_t{}.tif".format(input_dem_name.rstrip(".tif"), test_chunksize, workr, thrds)
        out_comb_vat_flat_path = os.path.abspath(os.path.join(output_dir_path, out_comb_vat_flat_name))

    compute_save_VAT_combined(general_combination, flat_combination, general_default, flat_default,
                              input_dem_path, out_comb_vat_path,
                             general_opacity, save_float, save_8bit, save_VAT_general,
                              out_comb_vat_general_path, save_VAT_flat, out_comb_vat_flat_path)


def compute_save_VAT_combined(general_combination, flat_combination, general_default, flat_default,
                              input_dem_path, out_comb_vat_path,
                              general_transparency, save_float, save_8bit, save_VAT_general, out_comb_vat_general_path,
                              save_VAT_flat, out_comb_vat_flat_path):
    dict_arr_res_nd = rvt.default.get_raster_dask_arr(raster_path=input_dem_path)  ## this is lazy load, option to:
   
    ## Option 1: compute entire array
    da_arr = dict_arr_res_nd['array']                                
    ## Option 2: or compute only desired region of the array, e.g. part of very large raster
    # h, w = dict_arr_res_nd['array'].shape
    # da_arr = dict_arr_res_nd['array'][0 : h//4, w//4: 3 * w//4]           

    # create and blend VAT general
    general_combination.add_dem_arr(dem_arr = da_arr,
                                    dem_resolution=dict_arr_res_nd["resolution"][0])
    if save_VAT_general:
        general_combination.add_dem_path(dem_path=input_dem_path)
        general_combination.render_all_images(save_render_path=out_comb_vat_general_path,
                                              save_float=save_float, save_8bit=save_8bit,
                                              default=general_default, no_data=dict_arr_res_nd["no_data"])
    # else:
    #     vat_arr_1 = general_combination.render_all_images(default=general_default, no_data=dict_arr_res_nd["no_data"])

    # create and blend VAT flat
    flat_combination.add_dem_arr(dem_arr = da_arr,
                                 dem_resolution=dict_arr_res_nd["resolution"][0])
    if save_VAT_flat:
        flat_combination.add_dem_path(dem_path=input_dem_path)
        flat_combination.render_all_images(save_render_path=out_comb_vat_flat_path,
                                           save_float=save_float, save_8bit=save_8bit,
                                           default=flat_default, no_data=dict_arr_res_nd["no_data"])
    # else:
    #     vat_arr_2 = flat_combination.render_all_images(default=flat_default, no_data=dict_arr_res_nd["no_data"])


    # # # create combination which blends VAT flat and VAT general together (check if flat or general are already computed and read from file)
    # if save_VAT_general and save_VAT_flat:
    #     combination = rvt.blend.BlenderCombination()
    #     image_path_1 = out_comb_vat_general_path
    #     image_path_2 = out_comb_vat_flat_path
    #     combination.create_layer(vis_method="VAT general", image_path=image_path_1, normalization="Value", minimum=0,
    #                             maximum=1, blend_mode="Normal", opacity=general_transparency)
    #     combination.create_layer(vis_method="VAT flat", image_path=image_path_2, normalization="Value", minimum=0,
    #                             maximum=1, blend_mode="Normal", opacity=100)

    #     combination.add_dem_path(dem_path=input_dem_path)

    #     # VAT combined
    #     combination.render_all_images(save_render_path=out_comb_vat_path, save_visualizations=False,
    #                                 save_float=save_float, save_8bit=save_8bit,
    #                                 no_data=dict_arr_res_nd["no_data"])
    #     out_comb_vat_8bit_path = out_comb_vat_path.rstrip("tif") + "_8bit.tif"
    #     if save_float and save_8bit:
    #         return "{} and {} successfully calculated and saved!".format(out_comb_vat_path, out_comb_vat_8bit_path)
    #     elif save_float:
    #         return "{} successfully calculated and saved!".format(out_comb_vat_path)
    #     elif save_8bit:
    #         return "{} successfully calculated and saved!".format(out_comb_vat_8bit_path)

def txt_output(output_file_path, output_file):
    with open(output_file_path, "w") as outf:
        json.dump(output_file, outf)
    return 0

def get_cluster_resources():
    """ Get number of cores and the amount of total memory available."""
    cluster = LocalCluster()
    client = Client(cluster)
    total_nr_cores = sum(client.ncores().values())            ## get total number of cores, do this once
    client.close()
    client.shutdown()
    total_memory = psutil.virtual_memory().total >> 30          ## bytes, shift to GB (aprox, rounding error)
    print(f"LocalCluster: total available cores: {total_nr_cores}, total memory: {total_memory} GB")

    return total_nr_cores, total_memory

# Program start
if __name__ == '__main__':
    total_cores, total_memory = get_cluster_resources()

    # total_cores = 16
    workers_threads_pairs = [(4, total_cores//4), (2, total_cores//2), (1, total_cores)]
    # workers_threads_pairs = [(1, total_cores - 2), (2, total_cores//2 - 2), (4, total_cores//4 - 2)] ### don't use all cores, for example - 2 threads per worker
    
    out=[]
    for n_workrs, n_thrds in workers_threads_pairs:
        for chunksize in chunksize_mib:
            print(f"Starting LocalCluster with: {n_workrs} workers and {n_thrds} threads per worker.\nArray chunksize: {chunksize}.")
            cluster = LocalCluster(n_workers = n_workrs, threads_per_worker = n_thrds,
                                    # memory_limit = f'{tot_memory // n_workrs}GB')    ## set memory limit (pauses workers)
                                    memory_limit = None)
            # dask_memusage.install(cluster.scheduler, f"memusage_{n_workrs}_{chunksize}.csv")
            client = Client(cluster)
            client.dashboard_link       ## <- open in browser to see dask diagnostics dashboard live
            start = time.time()

            try:
                dask.config.set({"array.chunk-size": chunksize})
                # print(dask.config.get("array.chunk-size"))
                combined_VAT(input_dir_path=input_dir_path, output_dir_path=output_dir_path,
                        test_chunksize =chunksize, workr = n_workrs, thrds = n_thrds, 
                        general_opacity=general_opacity, vat_combination_json_path=vat_combination_json_path,
                        terrains_sett_json_path=terrains_sett_json_path, save_float=save_float,
                        save_8bit=save_8bit, save_VAT_general=save_VAT_general, save_VAT_flat=save_VAT_flat)
                end = time.time()

                record = {'name': f"n_workers: {n_workrs}, n_threads: {n_thrds}", 
                        'duration': end - start, 
                        'chunk_size': chunksize, 
                        'n_workers': n_workrs, 
                        'n_threads': n_thrds}
                outfile = os.path.join(output_dir_path, f'w{n_workrs}_t{n_thrds}_{chunksize}.txt')  ###save intermediate results, in case of workers completly crashing 
                txt_output(output_file_path=outfile, output_file = record)

            except Exception as e:
                end = time.time()
                logger.exception(f"Exception {e}.")
                record = {'name': f"n_workers: {n_workrs}, n_threads: {n_thrds}", 
                        'exception at ': end - start,     ## save intermediate results, in case of workers completly crashing
                        'chunk_size': chunksize, 
                        'n_workers': n_workrs, 
                        'n_threads': n_thrds}
                outfile = os.path.join(output_dir_path, f'w{n_workrs}_t{n_thrds}_{chunksize}.txt')  
                txt_output(output_file_path=outfile, output_file = record)
                continue

            finally:
                out.append(record)
                print("Closing client.")
                client.close()
                client.shutdown()
                time.sleep(2)
                gc.collect()   #garbage collection, force (not sure if this does anything)

    outfile = os.path.join(output_dir_path, 'all_results.txt')  ###save all results joined
    txt_output(output_file_path=outfile, output_file = out)
