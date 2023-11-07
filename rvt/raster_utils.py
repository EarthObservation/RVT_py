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
    save_internal_nodata_mask: bool = False,
    internal_nodata_mask_array: Optional[npt.NDArray[Any]] = None,
) -> None:
    if save_internal_nodata_mask and internal_nodata_mask_array is None:
        raise ValueError(
            "If `save_internal_nodata_mask`=True, `internal_nodata_mask_array` needs to be provided!"
        )
    if not save_internal_nodata_mask and internal_nodata_mask_array is not None:
        raise ValueError(
            "If `internal_nodata_mask_array`=False `internal_nodata_mask_array` should be set to None!"
        )

    dataset_writer: rasterio.io.DatasetWriter
    with rasterio.Env() if save_internal_nodata_mask else nullcontext:  # if writing nodata mask array as
        with rasterio.open(output_path, "w", **rasterio_profile) as dataset_writer:
            dataset_writer.write(raster_array)
            if save_internal_nodata_mask:
                dataset_writer.dataset_mask(internal_nodata_mask_array)
