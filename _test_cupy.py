import cupy as cp
import rvt.visualizations
import rasterio as rio
import time
import numpy as np


def slope_aspect_gpu(
    dem,
    resolution_x,
    resolution_y,
    vertical_exaggeration_factor=1,
    output_units="radian",
):
    """
    Procedure can return terrain slope and aspect in radian units (default) or in alternative units (if specified).
    Slope is defined as 0 for Hz plane and pi/2 for vertical plane.
    Aspect iz defined as geographic azimuth: clockwise increasing, 0 or 2pi for the North direction.
    Currently applied finite difference method.

    Parameters
    ----------
    dem : input dem 2D cupy array
    resolution_x : dem resolution in X direction
    resolution_y : DEM resolution in Y direction
    vertical_exaggeration_factor : vertical exaggeration factor (must be greater than 0)
    output_units : percent, degree, radians

    Returns
    -------
    {"slope": slope_out, "aspect": aspect_out} : dictionaries with 2D numpy arrays
    """

    if vertical_exaggeration_factor <= 0:
        raise Exception(
            "rvt.vis.slope_aspect: vertical_exaggeration_factor must be a positive number!"
        )

    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.slope_aspect: resolution must be a positive number!")

    dem = dem.astype(np.float32)
    if vertical_exaggeration_factor != 1:
        dem = dem * vertical_exaggeration_factor

    # add frame of 0 (additional row up bottom and column left right)
    dem = cp.pad(dem, pad_width=1, mode="constant", constant_values=0)

    # derivatives in X and Y direction
    dzdx = ((cp.roll(dem, 1, axis=1) - cp.roll(dem, -1, axis=1)) / 2) / resolution_x
    dzdy = ((cp.roll(dem, -1, axis=0) - cp.roll(dem, 1, axis=0)) / 2) / resolution_y
    tan_slope = cp.sqrt(dzdx**2 + dzdy**2)

    # Compute slope
    if output_units == "percent":
        slope_out = tan_slope * 100
    elif output_units == "degree":
        slope_out = cp.rad2deg(np.arctan(tan_slope))
    elif output_units == "radian":
        slope_out = cp.arctan(tan_slope)
    else:
        raise Exception("rvt.vis.calculate_slope: Wrong function input 'output_units'!")

    # compute Aspect
    # aspect identifies the down slope direction of the maximum rate of change in value from each cell to its neighbors:
    #     0
    # 270    90
    #    180
    dzdy[
        dzdy == 0
    ] = 10e-9  # important for numeric stability - where dzdy is zero, make tangens to really high value

    aspect_out = cp.arctan2(dzdx, dzdy)  # atan2 took care of the quadrants
    if output_units == "degree":
        aspect_out = cp.rad2deg(aspect_out)

    # remove the frame (padding)
    slope_out = slope_out[1:-1, 1:-1]
    aspect_out = aspect_out[1:-1, 1:-1]

    # edges to -1
    slope_out[:, 0] = -1
    slope_out[0, :] = -1
    slope_out[:, -1] = -1
    slope_out[-1, :] = -1
    aspect_out[:, 0] = -1
    aspect_out[0, :] = -1
    aspect_out[:, -1] = -1
    aspect_out[-1, :] = -1

    return {"slope": slope_out, "aspect": aspect_out}


input_dem_path = r"test_data\TM1_564_146.tif"
input_dem_dataset = rio.open(input_dem_path)
t = input_dem_dataset.transform
x_res = t[0]
y_res = -t[4]
input_dem_arr = input_dem_dataset.read()[0]
input_dem_dataset.close()

print("DEM_path: {}, shape: {}".format(input_dem_path, input_dem_arr.shape))

# cupy
dem_cupy = cp.array(input_dem_arr)
start_time = time.time()
gpu_slope = slope_aspect_gpu(
    dem=dem_cupy,
    resolution_x=x_res,
    resolution_y=y_res,
    vertical_exaggeration_factor=1,
    output_units="degree",
)
end_time = time.time()
print("Gpu slope computing time: {}".format(end_time - start_time))
del gpu_slope
del dem_cupy

# numpy
dem_numpy = input_dem_arr
start_time = time.time()
cpu_slope = rvt.vis.slope_aspect(
    dem=dem_numpy,
    resolution_x=x_res,
    resolution_y=y_res,
    vertical_exaggeration_factor=1,
    output_units="degree",
)
end_time = time.time()
print("Cpu slope computing time: {}".format(end_time - start_time))
del cpu_slope
del dem_numpy
