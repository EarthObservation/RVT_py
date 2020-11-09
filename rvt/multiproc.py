import gdal
import rvt.vis
import rvt.default
import numpy as np
import multiprocessing as mp
import time


class RasterVisBlock:
    """Class for saving raster visualization block and its position in original raster."""
    def __init__(self, array, x, y):
        self.array = array
        self.x = x
        self.y = y


def save_multiprocess_vis(dem_path, vis, default, save_float=True, save_8bit=False, vis_float_path=None,
                          vis_8bit_path=None, x_block_size=5000, y_block_size=5000, offset=None):
    """
    Function reads and splits dem raster to blocks and parallel compute selected visualization on each block.
    Function is also parallel saving each block into output visualization raster. This function is suitable for really
    big rasters, 10000x10000 or bigger.

    Parameters
    ----------
    dem_path : str
        Path to input dem file.
    vis : str
        Selected visualization you want to compute: "hillshade, slope gradient, multiple directions hillshade,
        simple local relief model, sky-view factor, anisotropic sky-view factor, openness - positive,
        openness - negative, sky illumination, local dominance"
    default : rvt.default.DefaultValues
        Class where are all visualization parameters stored.
    save_float : bool
        If True it saves float raster.
    vis_float_path : str
        Path to output float raster visualization file (where you want to save output, file needs .tif ending).
    save_8bit : bool
        If Ture it saves 8bit raster (byte scaled float 0-255).
    vis_8bit_path : str
        Path to output 8bit raster visualization file (where you want to save output, file needs .tif ending).
    x_block_size : int
        Size of a block in x direction.
    y_block_size : int
        Size of a block in y direction
    offset : int
        Additional number of pixels on each side of a block when computing visualization (to remove block edges).
        If None function calculates it from default parameters.

    Returns
    -------
    Function saves output visualization GTiff in vis_path.
    """
    data_set = gdal.Open(dem_path)  # Open dem raster
    if data_set.RasterCount != 1:
        raise Exception("rvt.multiproc.save_multiprocess_vis: Input raster has more bands than 1!")
    gt = data_set.GetGeoTransform()
    x_res = gt[1]  # x_resolution
    y_res = -gt[5]  # y_resolution
    band = data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    del band

    nr_bands_float = 1
    nr_bands_8bit = 1
    if vis == "multiple directions hillshade":
        nr_bands_float = default.mhs_nr_dir
        nr_bands_8bit = 3

    if offset is None:  # if offset not set
        offset = get_block_offset(default=default, vis=vis)

    # create blank out raster
    if save_float:
        create_blank_raster(in_data_set=data_set, out_raster_path=vis_float_path, nr_bands=nr_bands_float, e_type=6)
    if save_8bit:
        create_blank_raster(in_data_set=data_set, out_raster_path=vis_8bit_path, nr_bands=nr_bands_8bit, e_type=1)

    # manager allows processes to all access same blocks_to_save list
    manager = mp.Manager()
    blocks_to_save_float = manager.list()  # list where all processes are saving computed float visualization blocks
    blocks_to_save_8bit = manager.list()  # list where all processes are saving computed 8bit visualization blocks

    # parallel process which is saving visualization blocks in output raster
    if save_float:
        save_float_process = mp.Process(target=save_block_visualization, args=(vis_float_path, blocks_to_save_float,
                                                                               nr_bands_float))
        save_float_process.start()
    if save_8bit:
        save_8bit_process = mp.Process(target=save_block_visualization, args=(vis_8bit_path, blocks_to_save_8bit,
                                                                               nr_bands_8bit))
        save_8bit_process.start()

    # list with processes
    blocks_processes = []

    # slice raster to blocks
    blocks = 0
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

            # create process for each block
            block_process = mp.Process(target=calc_add_vis_block, args=(save_float, save_8bit, blocks_to_save_float,
                                                                        blocks_to_save_8bit, block_array, default, x, y,
                                                                        left_offset, right_offset, top_offset,
                                                                        bottom_offset, vis, x_res, y_res))
            # start each process
            block_process.start()
            # add block to blocks_processes from where save_process is saving them
            blocks_processes.append(block_process)
            blocks += 1

    for block_process in blocks_processes:
        block_process.join()  # wait block_process to finish
    if save_float:
        blocks_to_save_float.append("kill")  # when block_process finished stop saving process
        save_float_process.join()
    if save_8bit:
        blocks_to_save_8bit.append("kill")  # when block_process finished stop saving process
        save_8bit_process.join()
    data_set = None  # close dem data set


def get_block_offset(default, vis):
    """Returns offset for selected visualization (vis). """
    if vis.lower() == "hillshade":
        return 1
    elif vis.lower() == "slope gradient":
        return 1
    elif vis.lower() == "multiple directions hillshade":
        return 1
    elif vis.lower() == "simple local relief model":
        return int(default.slrm_rad_cell)
    elif vis.lower() == "sky-view factor":
        return int(default.svf_r_max / 2)
    elif vis.lower() == "anisotropic sky-view factor":
        return int(default.svf_r_max / 2)
    elif vis.lower() == "openness - positive":
        return int(default.svf_r_max / 2)
    elif vis.lower() == "openness - negative":
        return int(default.svf_r_max / 2)
    elif vis.lower() == "sky illumination":
        return 10
    elif vis.lower() == "local dominance":
        return int(default.ld_max_rad)
    else:
        raise Exception("rvt.multiproc.calculate_block_visualization: Wrong vis parameter, vis options:"
                        "hillshade, slope gradient, multiple directions hillshade, simple local relief model,"
                        "sky-view factor, anisotropic sky-view factor, openness - positive, openness - negative, "
                        "sky illumination, local dominance")


def calc_add_vis_block(save_float, save_8bit, blocks_to_save_float, blocks_to_save_8bit, block_array, default, x, y,
                       left_offset, right_offset, top_offset, bottom_offset, vis, x_res, y_res):
    """Function calculates visualization on block, removes offset and adds block to list from where block is saved."""
    # calculate block visualization
    vis_block_array = calculate_block_visualization(dem_block_arr=block_array, vis=vis, default=default, x_res=x_res,
                                                    y_res=y_res)
    if vis != "multiple directions hillshade":
        del block_array

    # remove offset from visualization block
    if vis_block_array.ndim == 2:  # 2D array
        if right_offset == 0 and bottom_offset == 0:
            vis_block_array = vis_block_array[top_offset:, left_offset:]
        elif right_offset == 0:
            vis_block_array = vis_block_array[top_offset:-bottom_offset, left_offset:]
        elif bottom_offset == 0:
            vis_block_array = vis_block_array[top_offset:, left_offset :-right_offset]
        else:
            vis_block_array = vis_block_array[top_offset:-bottom_offset, left_offset:-right_offset]
    else:  # 3D array, multi hillshade
        if right_offset == 0 and bottom_offset == 0:
            vis_block_array = vis_block_array[:, top_offset:, left_offset:]
        elif right_offset == 0:
            vis_block_array = vis_block_array[:, top_offset:-bottom_offset, left_offset:]
        elif bottom_offset == 0:
            vis_block_array = vis_block_array[:, top_offset:, left_offset:-right_offset]
        else:
            vis_block_array = vis_block_array[:, top_offset:-bottom_offset, left_offset:-right_offset]

    if save_float:
        vis_block_float = RasterVisBlock(vis_block_array, x, y)
        blocks_to_save_float.append(vis_block_float)
    if save_8bit:
        if vis == "multiple directions hillshade":  # 8bit multihillshade (calculating hillshade in 3 directions)
            vis_block_array_8bit = default.float_to_8bit(float_arr=block_array, vis=vis, x_res=x_res, y_res=y_res)
            del block_array
            if right_offset == 0 and bottom_offset == 0:
                vis_block_array_8bit = vis_block_array_8bit[:, top_offset:, left_offset:]
            elif right_offset == 0:
                vis_block_array_8bit = vis_block_array_8bit[:, top_offset:-bottom_offset, left_offset:]
            elif bottom_offset == 0:
                vis_block_array_8bit = vis_block_array_8bit[:, top_offset:, left_offset:-right_offset]
            else:
                vis_block_array_8bit = vis_block_array_8bit[:, top_offset:-bottom_offset, left_offset:-right_offset]
        else:
            vis_block_array_8bit = default.float_to_8bit(float_arr=vis_block_array, vis=vis)
        vis_block_8bit = RasterVisBlock(vis_block_array_8bit, x, y)
        blocks_to_save_8bit.append(vis_block_8bit)

    del vis_block_array
    if save_8bit:
        del vis_block_array_8bit


def calculate_block_visualization(dem_block_arr, vis, default, x_res, y_res):
    """Calculate visualization on Block."""
    vis_arr = None
    if vis.lower() == "hillshade":
        vis_arr = default.get_hillshade(dem_arr=dem_block_arr, resolution_x=x_res, resolution_y=y_res)
    elif vis.lower() == "slope gradient":
        vis_arr = default.get_slope(dem_arr=dem_block_arr, resolution_x=x_res, resolution_y=y_res)["slope"]
    elif vis.lower() == "multiple directions hillshade":
        vis_arr = default.get_multi_hillshade(dem_arr=dem_block_arr, resolution_x=x_res, resolution_y=y_res)
    elif vis.lower() == "simple local relief model":
        vis_arr = default.get_slrm(dem_arr=dem_block_arr)
    elif vis.lower() == "sky-view factor":
        vis_arr = default.get_sky_view_factor(dem_arr=dem_block_arr, resolution=x_res, compute_svf=True,
                                              compute_asvf=False, compute_opns=False)["svf"]
    elif vis.lower() == "anisotropic sky-view factor":
        vis_arr = default.get_sky_view_factor(dem_arr=dem_block_arr, resolution=x_res, compute_svf=False,
                                              compute_asvf=True, compute_opns=False)["asvf"]
    elif vis.lower() == "openness - positive":
        vis_arr = default.get_sky_view_factor(dem_arr=dem_block_arr, resolution=x_res, compute_svf=False,
                                              compute_asvf=False, compute_opns=True)["opns"]
    elif vis.lower() == "openness - negative":
        vis_arr = default.get_neg_opns(dem_arr=dem_block_arr, resolution=x_res)
    elif vis.lower() == "sky illumination":
        vis_arr = default.get_sky_illumination(dem_arr=dem_block_arr, resolution=x_res)
    elif vis.lower() == "local dominance":
        vis_arr = default.get_local_dominance(dem_arr=dem_block_arr)
    else:
        raise Exception("rvt.multiproc.calculate_block_visualization: Wrong vis parameter, vis options:"
                        "hillshade, slope gradient, multiple directions hillshade, simple local relief model,"
                        "sky-view factor, anisotropic sky-view factor, openness - positive, openness - negative, "
                        "sky illumination, local dominance")
    return vis_arr


def save_block_visualization(vis_path, blocks_to_save, nr_bands=1):
    """Takes first block from blocks_to_save list and saves it to output raster then
     removes it form list."""
    out_data_set = gdal.Open(vis_path, gdal.GA_Update)
    while 1:
        if len(blocks_to_save) > 0:
            if blocks_to_save[0] == "kill":  # stops function/process
                out_data_set = None
                break
            if nr_bands == 1:  # one band
                out_data_set.GetRasterBand(1).WriteArray(blocks_to_save[0].array, blocks_to_save[0].x,
                                                         blocks_to_save[0].y)
                out_data_set.FlushCache()
            else:  # multiple bands
                for i_band in range(nr_bands):
                    band = i_band + 1
                    out_data_set.GetRasterBand(band).WriteArray(blocks_to_save[0].array[i_band], blocks_to_save[0].x,
                                                                blocks_to_save[0].y)
                    out_data_set.FlushCache()
            del blocks_to_save[0]  # remove block from blocks_to_save list


def create_blank_raster(in_data_set, out_raster_path, nr_bands=1, e_type=6):
    """Takes input data set and creates new raster. It copies input data set size, projection and geo info."""
    gtiff_driver = gdal.GetDriverByName("GTiff")
    band = in_data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    out_data_set = gtiff_driver.Create(out_raster_path, xsize=x_size, ysize=y_size, bands=nr_bands, eType=e_type,
                                       options=["BIGTIFF=IF_NEEDED"])
    out_data_set.SetProjection(in_data_set.GetProjection())
    out_data_set.SetGeoTransform(in_data_set.GetGeoTransform())
    out_data_set.FlushCache()
    out_data_set = None


# if __name__ == "__main__":
#     default = rvt.default.DefaultValues()
#     start_time = time.time()
#     x_block_size = 256
#     y_block_size = 256
#     save_multiprocess_vis(dem_path=r"..\test_data\TM1_564_146.tif",
#                           vis_float_path=r"..\test_data\TM1_564_146_multi_proc_tst.tif".format(
#                               x_block_size, y_block_size),
#                           vis_8bit_path=r"..\test_data\TM1_564_146_multi_proc_8bit_tst.tif".format(
#                               x_block_size, y_block_size),
#                           vis="local dominance",
#                           default=default,
#                           save_float=True,
#                           save_8bit=True,
#                           x_block_size=x_block_size,
#                           y_block_size=y_block_size)
#     end_time = time.time()
#     print("Multiprocessing: calculate and save hillshade with blocks({}x{}), dem(1000x1000), time={}".format(
#         x_block_size, y_block_size, end_time - start_time))

