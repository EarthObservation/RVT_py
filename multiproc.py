import gdal
import rvt.vis
import numpy as np
import multiprocessing as mp
import time


def save_multiprocess_vis(dem_path, vis_path, vis, offset, x_block_size, y_block_size):
    data_set = gdal.Open(dem_path)
    if data_set.RasterCount != 1:
        raise Exception("rvt.multiproc. : Input raster has more bands than 1!")
    gt = data_set.GetGeoTransform()
    x_res = gt[1]  # x_resolution
    y_res = -gt[5]  # y_resolution
    band = data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    del band
    blocks = 0

    # create out raster
    create_blank_raster(in_data_set=data_set, out_raster_path=vis_path)

    # multiprocess
    manager = mp.Manager()
    queue = manager.Queue()  # parameters: vis_block_arr, out_data_set, x, y
    pool = mp.Pool(mp.cpu_count() + 2)

    # listener for saving data
    watcher = pool.apply_async(save_block_visualization, (queue,))

    blocks_pools = []

    # slice raster to blocks
    for y in range(0, y_size, y_block_size):
        if y + y_block_size < y_size:  # if rows overlap
            rows = y_block_size
        else:
            rows = y_size - y
        for x in range(0, x_size, x_block_size):
            if x + x_block_size < x_size:  # if cols overlap
                cols = x_block_size
            else:
                cols = x_size - x

            # get offset for each block, check edges
            left_offset = 0
            right_offset = 0
            top_offset = 0
            bottom_offset = 0
            if x != 0:  # left offset
                if x - offset < 0:  # left overlap
                    left_offset = x
                else:
                    left_offset = offset
            if x + cols != x_size:  # right offset
                if x + cols + offset > x_size:  # right overlap
                    right_offset = x_size - x - cols
                else:
                    right_offset = offset
            if y != 0:  # apply top offset
                if y - offset < 0:  # top overlap
                    top_offset = y
                else:
                    top_offset = offset
            if y + rows != y_size:  # bottom offset
                if y + rows + offset > y_size:  # bottom overlap
                    bottom_offset = y_size - y - rows
                else:
                    bottom_offset = offset

            # reads block
            x_off = x - left_offset
            y_off = y - top_offset
            cols_off = cols + left_offset + right_offset
            rows_off = rows + top_offset + bottom_offset
            block_array = np.array(data_set.GetRasterBand(1).ReadAsArray(x_off, y_off, cols_off, rows_off))

            block_pool = pool.apply_async(calculate_and_save_block_vis, (block_array, vis_path, x, y,
                                                                         left_offset, right_offset, top_offset,
                                                                         bottom_offset, vis, x_res, y_res, queue))
            blocks_pools.append(block_pool)

            blocks += 1

    # collect results from the blocks through the pool result queue
    for block_pool in blocks_pools:
        block_pool.get()

    queue.put("kill")
    pool.close()
    pool.join()
    del data_set


def calculate_and_save_block_vis(block_array, vis_path, x, y, left_offset, right_offset, top_offset,
                                 bottom_offset, vis, x_res, y_res, queue):
    # calculate block visualization
    vis_block_array = calculate_block_visualization(dem_block_arr=block_array, vis=vis, x_res=x_res,
                                                    y_res=y_res)
    # remove offset from visualization block
    if right_offset == 0 and bottom_offset == 0:
        vis_block_array = vis_block_array[top_offset:, left_offset:]
    elif right_offset == 0:
        vis_block_array = vis_block_array[top_offset:-bottom_offset, left_offset:]
    elif bottom_offset == 0:
        vis_block_array = vis_block_array[top_offset:, left_offset:-right_offset]
    else:
        vis_block_array = vis_block_array[top_offset:-bottom_offset, left_offset:-right_offset]

    # save block
    # save_block_visualization(vis_block_arr=vis_block_array, out_data_set=out_data_set, x=x, y=y)
    queue.put([vis_block_array, vis_path, x, y])

    del vis_block_array


def calculate_block_visualization(dem_block_arr, vis, x_res, y_res):
    vis_arr = None
    if vis == "hillshade":
        vis_arr = rvt.vis.hillshade(dem=dem_block_arr, resolution_x=x_res, resolution_y=y_res)
    return vis_arr


def save_block_visualization(queue):
    # vis_block_arr, out_data_set, x, y
    while 1:
        list_of_arg = queue.get()
        # if kill stop saving data
        if list_of_arg == 'kill':
            break
        # read parameters from queue
        vis_block_arr = list_of_arg[0]
        vis_path = list_of_arg[1]
        x = list_of_arg[2]
        y = list_of_arg[3]
        try:
            out_data_set = gdal.Open(vis_path, gdal.GA_Update)
            out_data_set.GetRasterBand(1).WriteArray(vis_block_arr, x, y)
            out_data_set.FlushCache()
        except:  # TODO: solve it is not pythonic!
            time.sleep(0.2)
            out_data_set = gdal.Open(vis_path, gdal.GA_Update)
            out_data_set.GetRasterBand(1).WriteArray(vis_block_arr, x, y)
            out_data_set.FlushCache()
        del vis_block_arr
        del out_data_set
        del list_of_arg


def create_blank_raster(in_data_set, out_raster_path, bands=1, eType=6):
    """Takes input data set and creates new raster. It copies input data set size, projection and geo info."""
    gtiff_driver = gdal.GetDriverByName("GTiff")
    band = in_data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    out_data_set = gtiff_driver.Create(out_raster_path, xsize=x_size, ysize=y_size, bands=1, eType=6,
                                       options=["BIGTIFF=IF_NEEDED"])
    out_data_set.SetProjection(in_data_set.GetProjection())
    out_data_set.SetGeoTransform(in_data_set.GetGeoTransform())
    out_data_set.FlushCache()
    del out_data_set


if __name__ == "__main__":
    start_time = time.time()
    x_block_size = 1000
    y_block_size = 1000
    save_multiprocess_vis(dem_path=r"D:\RVT_py\test\multiproces_test\srtm_39and40_03.tif",
                          vis_path=r"D:\RVT_py\test\multiproces_test\srtm_39and40_03_mp{}x{}.tif".format(x_block_size,
                                                                                                         y_block_size),
                          vis="hillshade",
                          offset=5,
                          x_block_size=x_block_size,
                          y_block_size=y_block_size)
    end_time = time.time()
    print("Multiprocessing: calculate and save hillshade with blocks({}x{}), dem(12000x6000), time={}".format(
        x_block_size, y_block_size, end_time - start_time))
