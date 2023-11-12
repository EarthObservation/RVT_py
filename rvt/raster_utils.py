"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from contextlib import nullcontext
from pathlib import Path
from typing import Optional, Any
import rasterio
import rasterio.io
import rasterio.profiles
import numpy.typing as npt


def save_raster(
    output_path: Path,
    rasterio_profile: rasterio.profiles.Profile,
    raster_array: npt.NDArray[Any],
    internal_nodata_mask_array: Optional[npt.NDArray[Any]] = None,
) -> None:
    dataset_writer: rasterio.io.DatasetWriter
    with rasterio.Env() if internal_nodata_mask_array is not None else nullcontext:
        with rasterio.open(output_path, "w", **rasterio_profile) as dataset_writer:
            dataset_writer.write(raster_array)
            if internal_nodata_mask_array is not None:
                dataset_writer.dataset_mask(internal_nodata_mask_array)
