import gdal
import rvt.vis
import numpy as np


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
    band = None
    blocks = 0
    # create out raster
    gtiff_driver = gdal.GetDriverByName("GTiff")
    out_data_set = gtiff_driver.Create(vis_path, xsize=x_size, ysize=y_size, bands=1, eType=6,
                                       options=["BIGTIFF=IF_NEEDED"])
    out_data_set.SetProjection(data_set.GetProjection())
    out_data_set.SetGeoTransform(data_set.GetGeoTransform())

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

            calculate_and_save_block_vis(in_data_set=data_set, out_data_set=out_data_set, x=x, y=y, cols=cols,
                                         rows=rows, left_offset=left_offset, right_offset=right_offset,
                                         top_offset=top_offset, bottom_offset=bottom_offset,
                                         vis=vis, x_res=x_res, y_res=y_res)

            blocks += 1

    del data_set
    del out_data_set


def calculate_and_save_block_vis(in_data_set, out_data_set, x, y, cols, rows, left_offset, right_offset, top_offset,
                                 bottom_offset, vis, x_res, y_res):
    # reads block
    x_off = x - left_offset
    y_off = y - top_offset
    cols_off = cols + left_offset + right_offset
    rows_off = rows + top_offset + bottom_offset
    block_array = np.array(in_data_set.GetRasterBand(1).ReadAsArray(x_off, y_off, cols_off, rows_off))

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
    save_block_visualization(vis_block_arr=vis_block_array, out_data_set=out_data_set, x=x, y=y)

    del vis_block_array


def calculate_block_visualization(dem_block_arr, vis, x_res, y_res):
    vis_arr = None
    if vis == "hillshade":
        vis_arr = rvt.vis.hillshade(dem=dem_block_arr, resolution_x=x_res, resolution_y=y_res)
    return vis_arr


def save_block_visualization(vis_block_arr, out_data_set, x, y):
    out_data_set.GetRasterBand(1).WriteArray(vis_block_arr, x, y)
    out_data_set.FlushCache()
