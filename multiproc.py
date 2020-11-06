import gdal
import rvt.vis
import numpy as np
import multiprocessing as mp
import time


class rasterVisBlock:
    def __init__(self, array, x, y):
        self.array = array
        self.x = x
        self.y = y


def save_multiprocess_vis(dem_path, vis_path, vis, overlap, x_block_size, y_block_size):
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

    # listener for saving data
    manager = mp.Manager()
    blocks_to_save = manager.list()

    save_process = mp.Process(target=save_block_visualization, args=(vis_path, blocks_to_save))
    save_process.start()

    # list with processes
    blocks_processes = []

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

            # get overlap for each block, check edges
            left_offset = 0
            right_offset = 0
            top_offset = 0
            bottom_offset = 0
            if x != 0:  # left overlap
                if x - overlap < 0:  # left overlap
                    left_offset = x
                else:
                    left_offset = overlap
            if x + cols != x_size:  # right overlap
                if x + cols + overlap > x_size:  # right overlap
                    right_offset = x_size - x - cols
                else:
                    right_offset = overlap
            if y != 0:  # apply top overlap
                if y - overlap < 0:  # top overlap
                    top_offset = y
                else:
                    top_offset = overlap
            if y + rows != y_size:  # bottom overlap
                if y + rows + overlap > y_size:  # bottom overlap
                    bottom_offset = y_size - y - rows
                else:
                    bottom_offset = overlap

            # reads block
            x_off = x - left_offset
            y_off = y - top_offset
            cols_off = cols + left_offset + right_offset
            rows_off = rows + top_offset + bottom_offset
            block_array = np.array(data_set.GetRasterBand(1).ReadAsArray(x_off, y_off, cols_off, rows_off))

            block_process = mp.Process(target=calculate_and_save_block_vis, args=(blocks_to_save, block_array,
                                                                                  x, y, left_offset, right_offset,
                                                                                  top_offset, bottom_offset, vis, x_res,
                                                                                  y_res))
            block_process.start()
            blocks_processes.append(block_process)
            blocks += 1

    for block_process in blocks_processes:
        block_process.join()
    blocks_to_save.append("kill")
    save_process.join()
    data_set = None


def calculate_and_save_block_vis(blocks_to_save, block_array, x, y, left_offset, right_offset, top_offset,
                                 bottom_offset, vis, x_res, y_res):
    # calculate block visualization
    vis_block_array = calculate_block_visualization(dem_block_arr=block_array, vis=vis, x_res=x_res,
                                                    y_res=y_res)
    del block_array
    # remove overlap from visualization block
    if right_offset == 0 and bottom_offset == 0:
        vis_block_array = vis_block_array[top_offset:, left_offset:]
    elif right_offset == 0:
        vis_block_array = vis_block_array[top_offset:-bottom_offset, left_offset:]
    elif bottom_offset == 0:
        vis_block_array = vis_block_array[top_offset:, left_offset:-right_offset]
    else:
        vis_block_array = vis_block_array[top_offset:-bottom_offset, left_offset:-right_offset]

    vis_block = rasterVisBlock(vis_block_array, x, y)
    blocks_to_save.append(vis_block)

    del vis_block_array


def calculate_block_visualization(dem_block_arr, vis, x_res, y_res):
    vis_arr = None
    if vis == "hillshade":
        vis_arr = rvt.vis.hillshade(dem=dem_block_arr, resolution_x=x_res, resolution_y=y_res)
    return vis_arr


def save_block_visualization(vis_path, blocks_to_save):
    out_data_set = gdal.Open(vis_path, gdal.GA_Update)
    while 1:
        if len(blocks_to_save) > 0:
            if blocks_to_save[0] == "kill":
                out_data_set = None
                break
            # read parameters from queue
            out_data_set.GetRasterBand(1).WriteArray(blocks_to_save[0].array, blocks_to_save[0].x,
                                                     blocks_to_save[0].y)
            out_data_set.FlushCache()
            del blocks_to_save[0]


def create_blank_raster(in_data_set, out_raster_path, bands=1, e_type=6):
    """Takes input data set and creates new raster. It copies input data set size, projection and geo info."""
    gtiff_driver = gdal.GetDriverByName("GTiff")
    band = in_data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    out_data_set = gtiff_driver.Create(out_raster_path, xsize=x_size, ysize=y_size, bands=bands, eType=e_type,
                                       options=["BIGTIFF=IF_NEEDED"])
    out_data_set.SetProjection(in_data_set.GetProjection())
    out_data_set.SetGeoTransform(in_data_set.GetGeoTransform())
    out_data_set.FlushCache()
    out_data_set = None


if __name__ == "__main__":
    start_time = time.time()
    x_block_size = 2000
    y_block_size = 2000
    save_multiprocess_vis(dem_path=r"D:\RVT_py\test\multiproces_test\srtm_39and40_03.tif",
                          vis_path=r"D:\RVT_py\test\multiproces_test\srtm_39and40_03_mp{}x{}.tif".format(x_block_size,
                                                                                                         y_block_size),
                          vis="hillshade",
                          overlap=5,
                          x_block_size=x_block_size,
                          y_block_size=y_block_size)
    end_time = time.time()
    print("Multiprocessing: calculate and save hillshade with blocks({}x{}), dem(12000x6000), time={}".format(
        x_block_size, y_block_size, end_time - start_time))