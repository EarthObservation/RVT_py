"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass, fields
from enum import Enum
from pathlib import Path
from typing import Optional, Tuple, Any, Dict, Type

import numpy.typing as npt
import json
import numpy as np
import rasterio
import rasterio.io
import rasterio.profiles
import rasterio.enums
from overrides import override

from rvt.blender_functions import Normalization
from rvt.enums import (
    RVTVisualizationName,
    OpennessType,
    SlopeOutputUnit,
    SvfNoiseRemove,
    AnisotropyLevel,
    NormalizationMode,
)
from rvt.raster_utils import save_raster
from rvt.visualizations import (
    slope_aspect,
    hillshade,
    multi_hillshade,
    simple_local_relief_model,
    horizon_visualizations,
    HorizonVisualizationResult,
    local_dominance,
    multi_scale_relief_model,
    multi_scale_topographic_position,
    shadow_horizon,
    byte_scale,
)

_HORIZON_VISUALIZATIONS_DEFAULT_NUMBER_OF_DIRECTIONS: int = 16
_HORIZON_VISUALIZATIONS_DEFAULT_MAXIMUM_SEARCH_RADIUS: int = 10
_HORIZON_VISUALIZATIONS_DEFAULT_NOISE_REMOVE: SvfNoiseRemove = SvfNoiseRemove.NO_REMOVE
_HORIZON_VISUALIZATIONS_DEFAULT_DIRECTION_OF_ANISOTROPY: float = 315.0
_HORIZON_VISUALIZATIONS_DEFAULT_ANISOTROPY_LEVEL: AnisotropyLevel = AnisotropyLevel.LOW


class To8bit(Normalization):
    def to_string(self) -> str:
        return f"normalization_mode={self.normalization_mode.value}, min={self.minimum}, max={self.maximum}"

    @classmethod
    def from_string(cls, string: str) -> To8bit:
        elements = string.split(",")
        return To8bit(
            normalization_mode=NormalizationMode(
                elements[0].strip().lstrip("normalization_mode=")
            ),
            minimum=float(elements[1].strip().lstrip("min=")),
            maximum=float(elements[2].strip().lstrip("max=")),
        )

    def from_float(self, float_image: npt.NDArray[Any], no_data: Optional[float]):
        return byte_scale(
            data=self.normalize_image(image=float_image),
            c_min=0,
            c_max=1,
            high=255,
            low=0,
            no_data=no_data,
        )


class RVTVisualization(ABC):
    RVT_VISUALIZATION_NAME: RVTVisualizationName
    to_8bit: To8bit

    @abstractmethod
    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        """
        Get visualization file name with parameters from DEM file name. Function returns visualization name with
        tif suffix.
        """
        pass

    def get_visualization_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return directory_path / self.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    @abstractmethod
    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        pass

    def convert_visualization_to_8bit(
        self, visualization_array: npt.NDArray[Any], no_data: Optional[float] = np.nan
    ) -> npt.NDArray[Any]:
        if len(visualization_array.shape) == 3 and visualization_array.shape[2] != 3:
            raise RuntimeError(
                f"Can't convert visualization array to 8bit, invalid shape ({visualization_array.shape}). "
                "Only single band (2D array) or 3 band (3D array) supported."
            )
        return self.to_8bit.from_float(float_image=visualization_array, no_data=no_data)

    def compute_visualization_from_dem_file(
        self, dem_path: Path, vertical_exaggeration_factor: float
    ) -> npt.NDArray[Any]:
        dem_file: rasterio.io.DatasetReader
        with rasterio.open(dem_path, "r") as dem_file:
            return self.compute_visualization(
                dem=dem_file.read(0),
                dem_resolution=(abs(dem_file.transform[0]) + abs(dem_file.transform[4]))
                / 2,
                dem_nodata=dem_file.nodata,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
            )

    @classmethod
    def _validate_save_visualization_parameters(
        cls,
        dem_path: Path,
        overwrite: bool,
        save_float_visualization: bool,
        output_float_visualization_path: Optional[Path],
        save_8bit_visualization: bool,
        output_8bit_visualization_path: Optional[Path],
    ) -> None:
        if not dem_path.exists():
            raise ValueError(f"Dem path {dem_path} doesn't exist!")
        if not save_float_visualization and not save_8bit_visualization:
            raise ValueError(
                "In order to save visualization `save_float_visualization` or `save_8bit_visualization` "
                "needs to be True!"
            )
        if save_float_visualization and output_float_visualization_path is None:
            raise ValueError(
                "If `save_float_visualization`=True, `output_float_visualization_path` needs to be provided!"
            )
        if save_8bit_visualization and output_8bit_visualization_path is None:
            raise ValueError(
                "If `save_8bit_visualization`=True, `output_8bit_visualization_path` needs to be provided!"
            )
        if not overwrite:
            if save_float_visualization and output_float_visualization_path.exists():
                raise ValueError(
                    f"Output visualization path ({output_float_visualization_path}) already exists!"
                )
            if save_8bit_visualization and output_8bit_visualization_path.exists():
                raise ValueError(
                    f"Output visualization path ({output_8bit_visualization_path}) already exists!"
                )

    def save_visualization(
        self,
        dem_path: Path,
        vertical_exaggeration_factor: float,
        overwrite: bool = True,
        compression: rasterio.enums.Compression = rasterio.enums.Compression.lzw,
        save_float_visualization: bool = True,
        output_float_visualization_path: Optional[Path] = None,
        output_float_visualization_nodata: Optional[float] = np.nan,
        save_8bit_visualization: bool = False,
        output_8bit_visualization_path: Optional[Path] = None,
    ) -> None:
        self._validate_save_visualization_parameters(
            dem_path=dem_path,
            overwrite=overwrite,
            save_float_visualization=save_float_visualization,
            output_float_visualization_path=output_float_visualization_path,
            save_8bit_visualization=save_8bit_visualization,
            output_8bit_visualization_path=output_8bit_visualization_path,
        )
        dem_file: rasterio.io.DatasetReader
        with rasterio.open(dem_path, "r") as dem_file:
            dem_transform = dem_file.profile["transform"]
            visualization_array = self.compute_visualization(
                dem=dem_file.read(0),
                dem_resolution=(abs(dem_file.transform[0]) + abs(dem_file.transform[4]))
                / 2,
                dem_nodata=dem_file.nodata,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
            )
        number_of_bands = (
            1 if len(visualization_array.shape) == 2 else visualization_array.shape[0]
        )
        if save_float_visualization:
            if not np.isnan(output_float_visualization_nodata):  # change nodata value
                visualization_array[
                    np.isnan(visualization_array)
                ] = output_float_visualization_nodata
            visualization_profile = rasterio.profiles.DefaultGTiffProfile(
                transform=dem_transform,
                compression=compression.value,
                count=number_of_bands,
                dtype=rasterio.float32,
                nodata=output_float_visualization_nodata,
            )
            save_raster(
                output_path=output_float_visualization_path,
                rasterio_profile=visualization_profile,
                raster_array=visualization_array,
            )

        if save_8bit_visualization:
            visualization_8bit_array = self.convert_visualization_to_8bit(
                visualization_array=visualization_array, no_data=np.nan
            )
            nodata_mask = visualization_array[np.isnan(visualization_array)]
            visualization_8bit_array[nodata_mask] = 0
            visualization_profile = rasterio.profiles.DefaultGTiffProfile(
                transform=dem_transform,
                compression=compression.value,
                count=number_of_bands,
                dtype=rasterio.uint8,
                nodata=None,
            )
            save_raster(
                output_path=output_8bit_visualization_path,
                rasterio_profile=visualization_profile,
                raster_array=visualization_8bit_array,
                internal_nodata_mask_array=nodata_mask,
            )


class Slope(RVTVisualization):
    RVT_VISUALIZATION_NAME: RVTVisualizationName = RVTVisualizationName.SLOPE
    _DEFAULT_TO_8BIT: To8bit = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.00, maximum=51.00
    )
    _DEFAULT_OUTPUT_UNIT: SlopeOutputUnit = SlopeOutputUnit.DEGREE

    def __init__(
        self,
        output_unit: SlopeOutputUnit = _DEFAULT_OUTPUT_UNIT,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.output_unit = output_unit
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = f"{dem_file_name}_SLOPE"
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return slope_aspect(
            dem=dem,
            resolution_x=dem_resolution,
            resolution_y=dem_resolution,
            output_unit=self.output_unit,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        ).slope


class Shadow(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.SHADOW
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.00, maximum=1.00
    )
    _DEFAULT_SUN_AZIMUTH = 315.0
    _DEFAULT_SUN_ELEVATION = 35.0

    def __init__(
        self,
        sun_azimuth: float = _DEFAULT_SUN_AZIMUTH,
        sun_elevation: float = _DEFAULT_SUN_ELEVATION,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.sun_azimuth = sun_azimuth
        self.sun_elevation = sun_elevation
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_SHADOW_A{self.sun_azimuth}_H{self.sun_elevation}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return shadow_horizon(
            dem=dem,
            resolution=dem_resolution,
            shadow_az=self.sun_azimuth,
            shadow_el=self.sun_elevation,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )["shadow"]


class Hillshade(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.HILLSHADE
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.00, maximum=1.00
    )
    _DEFAULT_SUN_AZIMUTH = 315.0
    _DEFAULT_SUN_ELEVATION = 35.0

    def __init__(
        self,
        sun_azimuth: float = _DEFAULT_SUN_AZIMUTH,
        sun_elevation: float = _DEFAULT_SUN_ELEVATION,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.sun_azimuth = sun_azimuth
        self.sun_elevation = sun_elevation
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = f"{dem_file_name}_HS_A{self.sun_azimuth}_H{self.sun_elevation}"
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return hillshade(
            dem=dem,
            resolution_x=dem_resolution,
            resolution_y=dem_resolution,
            sun_azimuth=self.sun_azimuth,
            sun_elevation=self.sun_elevation,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )


class MultipleDirectionsHillshade(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.00, maximum=1.00
    )
    _DEFAULT_NUMBER_OF_DIRECTIONS: int = 16
    _DEFAULT_SUN_ELEVATION: float = 35.0
    _SUN_AZIMUTHS_FOR_8BIT_VISUALIZATION = (315, 22.5, 90)

    def __init__(
        self,
        number_of_directions: int = _DEFAULT_NUMBER_OF_DIRECTIONS,
        sun_elevation: float = _DEFAULT_SUN_ELEVATION,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.number_of_directions = number_of_directions
        self.sun_elevation = sun_elevation
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = f"{dem_file_name}_MULTI-HS_D{self.number_of_directions}_H{self.sun_elevation}"
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return multi_hillshade(
            dem=dem,
            resolution_x=dem_resolution,
            resolution_y=dem_resolution,
            number_of_directions=self.number_of_directions,
            sun_elevation=self.sun_elevation,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )

    @override
    def convert_visualization_to_8bit(
        self, visualization_array: npt.NDArray[Any], no_data: Optional[float] = np.nan
    ) -> npt.NDArray[Any]:
        raise NotImplementedError(
            "Impossible to convert Multiple directions hillshade float visualization to 8bit. "
            "Use `compute_8bit_visualization()` method instead."
        )

    def compute_8bit_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        red_channel = self.to_8bit.from_float(
            float_image=Hillshade(
                sun_azimuth=self._SUN_AZIMUTHS_FOR_8BIT_VISUALIZATION[0],
                sun_elevation=self.sun_elevation,
            ).compute_visualization(
                dem=dem,
                dem_resolution=dem_resolution,
                dem_nodata=dem_nodata,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
            ),
            no_data=np.nan,
        )
        green_channel = self.to_8bit.from_float(
            float_image=Hillshade(
                sun_azimuth=self._SUN_AZIMUTHS_FOR_8BIT_VISUALIZATION[1],
                sun_elevation=self.sun_elevation,
            ).compute_visualization(
                dem=dem,
                dem_resolution=dem_resolution,
                dem_nodata=dem_nodata,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
            ),
            no_data=np.nan,
        )
        blue_channel = self.to_8bit.from_float(
            float_image=Hillshade(
                sun_azimuth=self._SUN_AZIMUTHS_FOR_8BIT_VISUALIZATION[2],
                sun_elevation=self.sun_elevation,
            ).compute_visualization(
                dem=dem,
                dem_resolution=dem_resolution,
                dem_nodata=dem_nodata,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
            ),
            no_data=np.nan,
        )
        return np.array([red_channel, green_channel, blue_channel])

    @override
    def save_visualization(
        self,
        dem_path: Path,
        vertical_exaggeration_factor: float,
        overwrite: bool = True,
        compression: rasterio.enums.Compression = rasterio.enums.Compression.lzw,
        save_float_visualization: bool = True,
        output_float_visualization_path: Optional[Path] = None,
        output_float_visualization_nodata: Optional[float] = np.nan,
        save_8bit_visualization: bool = False,
        output_8bit_visualization_path: Optional[Path] = None,
    ) -> None:
        if not save_float_visualization and not save_8bit_visualization:
            raise ValueError(
                "In order to save visualization `save_float_visualization` or `save_8bit_visualization` "
                "needs to be True!"
            )
        if save_8bit_visualization and output_8bit_visualization_path is None:
            raise ValueError(
                "If `save_8bit_visualization`=True, `output_8bit_visualization_path` needs to be provided!"
            )
        if not overwrite:
            if save_8bit_visualization and output_8bit_visualization_path.exists():
                raise ValueError(
                    f"Output visualization path ({output_8bit_visualization_path}) already exists!"
                )
        if save_float_visualization:
            super().save_visualization(
                dem_path=dem_path,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
                overwrite=overwrite,
                compression=compression,
                save_float_visualization=True,
                output_float_visualization_path=output_float_visualization_path,
                output_float_visualization_nodata=output_float_visualization_nodata,
                save_8bit_visualization=False,
                output_8bit_visualization_path=None,
            )
        if save_8bit_visualization:
            with rasterio.open(dem_path) as dem_file:
                visualization_8bit_array = self.compute_8bit_visualization(
                    dem=dem_file.read(0),
                    dem_resolution=(
                        abs(dem_file.transform[0]) + abs(dem_file.transform[4])
                    )
                    / 2,
                    dem_nodata=dem_file.nodata,
                    vertical_exaggeration_factor=vertical_exaggeration_factor,
                )
                nodata_mask = visualization_8bit_array[
                    np.isnan(visualization_8bit_array)
                ]
                visualization_8bit_array[nodata_mask] = 0
                visualization_profile = rasterio.profiles.DefaultGTiffProfile(
                    transform=dem_file.profile["transform"],
                    compression=compression.value,
                    count=3,
                    dtype=rasterio.uint8,
                    nodata=None,
                )
                save_raster(
                    output_path=output_8bit_visualization_path,
                    rasterio_profile=visualization_profile,
                    raster_array=visualization_8bit_array,
                    internal_nodata_mask_array=nodata_mask,
                )


class SimpleLocalReliefModel(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.SIMPLE_LOCAL_RELIEF_MODEL
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=-2.00, maximum=2.00
    )
    _DEFAULT_RADIUS_CELL: int = 20

    def __init__(
        self,
        radius_cell: int = _DEFAULT_RADIUS_CELL,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.radius_cell = radius_cell
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = f"{dem_file_name}_SLRM_R{self.radius_cell}"
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return simple_local_relief_model(
            dem=dem,
            radius_cell=self.radius_cell,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )


class SkyViewFactor(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.SKY_VIEW_FACTOR
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.6375, maximum=1.00
    )
    _DEFAULT_NUMBER_OF_DIRECTIONS: int = (
        _HORIZON_VISUALIZATIONS_DEFAULT_NUMBER_OF_DIRECTIONS
    )
    _DEFAULT_MAXIMUM_SEARCH_RADIUS: int = (
        _HORIZON_VISUALIZATIONS_DEFAULT_MAXIMUM_SEARCH_RADIUS
    )
    _DEFAULT_NOISE_REMOVE: SvfNoiseRemove = _HORIZON_VISUALIZATIONS_DEFAULT_NOISE_REMOVE

    def __init__(
        self,
        number_of_directions: int = _DEFAULT_NUMBER_OF_DIRECTIONS,
        maximum_search_radius: int = _DEFAULT_MAXIMUM_SEARCH_RADIUS,
        noise_remove: SvfNoiseRemove = _DEFAULT_NOISE_REMOVE,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.number_of_directions = number_of_directions
        self.maximum_search_radius = maximum_search_radius
        self.noise_remove = noise_remove
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_SVF_R{self.maximum_search_radius}_D{self.number_of_directions}_"
            f"NR{self.noise_remove.value}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return horizon_visualizations(
            dem=dem,
            resolution=dem_resolution,
            compute_svf=True,
            compute_opns=False,
            compute_asvf=False,
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        ).sky_view_factor


class AnisotropicSkyViewFactor(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.ANISOTROPIC_SKY_VIEW_FACTOR
    _DEFAULT_TO_8BIT: To8bit = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.70, maximum=0.90
    )
    _DEFAULT_NUMBER_OF_DIRECTIONS: int = (
        _HORIZON_VISUALIZATIONS_DEFAULT_NUMBER_OF_DIRECTIONS
    )
    _DEFAULT_MAXIMUM_SEARCH_RADIUS: int = (
        _HORIZON_VISUALIZATIONS_DEFAULT_MAXIMUM_SEARCH_RADIUS
    )
    _DEFAULT_NOISE_REMOVE: SvfNoiseRemove = _HORIZON_VISUALIZATIONS_DEFAULT_NOISE_REMOVE
    _DEFAULT_DIRECTION_OF_ANISOTROPY: float = (
        _HORIZON_VISUALIZATIONS_DEFAULT_DIRECTION_OF_ANISOTROPY
    )
    _DEFAULT_ANISOTROPY_LEVEL: AnisotropyLevel = (
        _HORIZON_VISUALIZATIONS_DEFAULT_ANISOTROPY_LEVEL
    )

    def __init__(
        self,
        number_of_directions: int = _DEFAULT_NUMBER_OF_DIRECTIONS,
        maximum_search_radius: int = _DEFAULT_MAXIMUM_SEARCH_RADIUS,
        noise_remove: SvfNoiseRemove = _DEFAULT_NOISE_REMOVE,
        direction_of_anisotropy: float = _DEFAULT_DIRECTION_OF_ANISOTROPY,
        anisotropy_level: AnisotropyLevel = _DEFAULT_ANISOTROPY_LEVEL,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.number_of_directions = number_of_directions
        self.maximum_search_radius = maximum_search_radius
        self.noise_remove = noise_remove
        self.direction_of_anisotropy = direction_of_anisotropy
        self.anisotropy_level = anisotropy_level
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_SVF-A_R{self.maximum_search_radius}_D{self.number_of_directions}_"
            f"NR{self.noise_remove.value}_A{self.direction_of_anisotropy}_AL{self.anisotropy_level.value}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return horizon_visualizations(
            dem=dem,
            resolution=dem_resolution,
            compute_svf=False,
            compute_opns=False,
            compute_asvf=True,
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
            direction_of_anisotropy=self.direction_of_anisotropy,
            anisotropy_level=self.anisotropy_level,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        ).anisotropic_sky_view_factor


class Openness(RVTVisualization):
    @property
    def RVT_VISUALIZATION_NAME(self) -> RVTVisualizationName:
        if self.openness_type == OpennessType.POSITIVE:
            return RVTVisualizationName.POSITIVE_OPENNESS
        elif self.openness_type == OpennessType.NEGATIVE:
            return RVTVisualizationName.NEGATIVE_OPENNESS
        else:
            raise ValueError(f"Wrong Openness type {self.openness_type}!")

    _DEFAULT_POSITIVE_OPENNESS_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=60.00, maximum=95.00
    )
    _DEFAULT_NEGATIVE_OPENNESS_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=60.00, maximum=95.00
    )
    _DEFAULT_OPENNESS_TYPE: OpennessType = OpennessType.POSITIVE
    _DEFAULT_NUMBER_OF_DIRECTIONS: int = (
        _HORIZON_VISUALIZATIONS_DEFAULT_NUMBER_OF_DIRECTIONS
    )
    _DEFAULT_MAXIMUM_SEARCH_RADIUS: int = (
        _HORIZON_VISUALIZATIONS_DEFAULT_MAXIMUM_SEARCH_RADIUS
    )
    _DEFAULT_NOISE_REMOVE: SvfNoiseRemove = _HORIZON_VISUALIZATIONS_DEFAULT_NOISE_REMOVE

    def __init__(
        self,
        openness_type: OpennessType = _DEFAULT_OPENNESS_TYPE,
        number_of_directions: int = _DEFAULT_NUMBER_OF_DIRECTIONS,
        maximum_search_radius: int = _DEFAULT_MAXIMUM_SEARCH_RADIUS,
        noise_remove: SvfNoiseRemove = _DEFAULT_NOISE_REMOVE,
        to_8bit: Optional[To8bit] = None,
    ):
        self.openness_type = openness_type
        self.number_of_directions = number_of_directions
        self.maximum_search_radius = maximum_search_radius
        self.noise_remove = noise_remove
        if to_8bit is None:
            if openness_type == OpennessType.POSITIVE:
                self.to_8bit = self._DEFAULT_POSITIVE_OPENNESS_TO_8BIT
            if openness_type == OpennessType.NEGATIVE:
                self.to_8bit = self._DEFAULT_NEGATIVE_OPENNESS_TO_8BIT
        else:
            self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_OPEN-{self.openness_type.name[:3]}_R{self.maximum_search_radius}_"
            f"D{self.number_of_directions}_NR{self.noise_remove.value}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        if OpennessType.NEGATIVE:
            dem = -1 * dem
        return horizon_visualizations(
            dem=dem,
            resolution=dem_resolution,
            compute_svf=True,
            compute_opns=False,
            compute_asvf=False,
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        ).openness


@dataclass
class HorizonVisualizations:
    number_of_directions: int = _HORIZON_VISUALIZATIONS_DEFAULT_NUMBER_OF_DIRECTIONS
    maximum_search_radius: int = _HORIZON_VISUALIZATIONS_DEFAULT_MAXIMUM_SEARCH_RADIUS
    noise_remove: SvfNoiseRemove = _HORIZON_VISUALIZATIONS_DEFAULT_NOISE_REMOVE
    direction_of_anisotropy: float = (
        _HORIZON_VISUALIZATIONS_DEFAULT_DIRECTION_OF_ANISOTROPY
    )
    anisotropy_level: AnisotropyLevel = _HORIZON_VISUALIZATIONS_DEFAULT_ANISOTROPY_LEVEL

    @property
    def sky_view_factor(self) -> SkyViewFactor:
        return SkyViewFactor(
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
        )

    @property
    def anisotropic_sky_view_factor(self) -> AnisotropicSkyViewFactor:
        return AnisotropicSkyViewFactor(
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
            direction_of_anisotropy=self.direction_of_anisotropy,
            anisotropy_level=self.anisotropy_level,
        )

    @property
    def positive_openness(self) -> Openness:
        return Openness(
            openness_type=OpennessType.POSITIVE,
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
        )

    @property
    def negative_openness(self) -> Openness:
        return Openness(
            openness_type=OpennessType.NEGATIVE,
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
        )

    def compute_horizon_visualizations(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
        compute_svf: bool,
        compute_opns: bool,
        compute_asvf: bool,
    ) -> HorizonVisualizationResult:
        return horizon_visualizations(
            dem=dem,
            resolution=dem_resolution,
            compute_svf=compute_svf,
            compute_opns=compute_opns,
            compute_asvf=compute_asvf,
            number_of_directions=self.number_of_directions,
            maximum_search_radius=self.maximum_search_radius,
            noise_remove=self.noise_remove,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )

    @classmethod
    def _validate_save_horizon_visualization_parameters(
        cls,
        dem_path: Path,
        overwrite: bool,
        save_svf_float_visualization: bool,
        output_svf_float_visualization_path: Optional[Path],
        save_svf_8bit_visualization: bool,
        output_svf_8bit_visualization_path: Optional[Path],
        save_asvf_float_visualization: bool,
        output_asvf_float_visualization_path: Optional[Path],
        save_asvf_8bit_visualization: bool,
        output_asvf_8bit_visualization_path: Optional[Path],
        save_opns_float_visualization: bool,
        output_opns_float_visualization_path: Optional[Path],
        save_opns_8bit_visualization: bool,
        output_opns_8bit_visualization_path: Optional[Path],
    ) -> None:
        if not dem_path.exists():
            raise ValueError(f"Dem path {dem_path} doesn't exist!")
        if not save_svf_float_visualization and not save_svf_8bit_visualization:
            raise ValueError(
                "In order to save visualization one of the `save_{vis}_float_visualization` or"
                " `save_{vis}_8bit_visualization` needs to be True!"
            )
        if save_svf_float_visualization and output_svf_float_visualization_path is None:
            raise ValueError(
                "If `save_svf_float_visualization`=True, `output_svf_float_visualization_path` needs to be provided!"
            )
        if save_svf_8bit_visualization and output_svf_8bit_visualization_path is None:
            raise ValueError(
                "If `save_svf_8bit_visualization`=True, `output_svf_8bit_visualization_path` needs to be provided!"
            )
        if (
            save_asvf_float_visualization
            and output_asvf_float_visualization_path is None
        ):
            raise ValueError(
                "If `save_asvf_float_visualization`=True, `output_asvf_float_visualization_path` needs to be provided!"
            )
        if save_asvf_8bit_visualization and output_asvf_8bit_visualization_path is None:
            raise ValueError(
                "If `save_asvf_8bit_visualization`=True, `output_asvf_8bit_visualization_path` needs to be provided!"
            )
        if (
            save_opns_float_visualization
            and output_opns_float_visualization_path is None
        ):
            raise ValueError(
                "If `save_opns_float_visualization`=True, `output_opns_float_visualization_path` needs to be provided!"
            )
        if save_opns_8bit_visualization and output_opns_8bit_visualization_path is None:
            raise ValueError(
                "If `save_opns_8bit_visualization`=True, `output_opns_8bit_visualization_path` needs to be provided!"
            )
        if not overwrite:
            if (
                save_svf_float_visualization
                and output_svf_float_visualization_path.exists()
            ):
                raise ValueError(
                    f"Output visualization path ({output_svf_float_visualization_path}) already exists!"
                )
            if (
                save_svf_8bit_visualization
                and output_svf_8bit_visualization_path.exists()
            ):
                raise ValueError(
                    f"Output visualization path ({output_svf_8bit_visualization_path}) already exists!"
                )
            if (
                save_asvf_float_visualization
                and output_asvf_float_visualization_path.exists()
            ):
                raise ValueError(
                    f"Output visualization path ({output_asvf_float_visualization_path}) already exists!"
                )
            if (
                save_asvf_8bit_visualization
                and output_asvf_8bit_visualization_path.exists()
            ):
                raise ValueError(
                    f"Output visualization path ({output_asvf_8bit_visualization_path}) already exists!"
                )
            if (
                save_opns_float_visualization
                and output_opns_float_visualization_path.exists()
            ):
                raise ValueError(
                    f"Output visualization path ({output_opns_float_visualization_path}) already exists!"
                )
            if (
                save_opns_8bit_visualization
                and output_opns_8bit_visualization_path.exists()
            ):
                raise ValueError(
                    f"Output visualization path ({output_opns_8bit_visualization_path}) already exists!"
                )

    def save_horizon_visualizations(
        self,
        dem_path: Path,
        vertical_exaggeration_factor: float,
        overwrite: bool = True,
        compression: rasterio.enums.Compression = rasterio.enums.Compression.lzw,
        save_svf_float_visualization: bool = True,
        output_svf_float_visualization_path: Optional[Path] = None,
        save_svf_8bit_visualization: bool = False,
        output_svf_8bit_visualization_path: Optional[Path] = None,
        save_asvf_float_visualization: bool = True,
        output_asvf_float_visualization_path: Optional[Path] = None,
        save_asvf_8bit_visualization: bool = False,
        output_asvf_8bit_visualization_path: Optional[Path] = None,
        save_opns_float_visualization: bool = True,
        output_opns_float_visualization_path: Optional[Path] = None,
        save_opns_8bit_visualization: bool = False,
        output_opns_8bit_visualization_path: Optional[Path] = None,
        output_float_visualization_nodata: Optional[float] = np.nan,
    ):
        def _save_float_visualization(
            visualization_array: npt.NDArray[Any], output_path: Path
        ) -> None:
            assert visualization_array is not None
            if not np.isnan(output_float_visualization_nodata):  # change nodata value
                visualization_array = visualization_array[
                    np.isnan(svf_visualization_array)
                ] = output_float_visualization_nodata
            number_of_bands = 1
            visualization_profile = rasterio.profiles.DefaultGTiffProfile(
                transform=dem_transform,
                compression=compression.value,
                count=number_of_bands,
                dtype=rasterio.float32,
                nodata=output_float_visualization_nodata,
            )
            save_raster(
                output_path=output_path,
                rasterio_profile=visualization_profile,
                raster_array=visualization_array,
            )

        def _save_8bit_visualization(
            visualization_8bit_array: npt.NDArray[Any], output_path: Path
        ) -> None:
            nodata_mask = visualization_8bit_array[np.isnan(visualization_8bit_array)]
            visualization_8bit_array[nodata_mask] = 0
            visualization_profile = rasterio.profiles.DefaultGTiffProfile(
                transform=dem_transform,
                compression=compression.value,
                count=1,
                dtype=rasterio.uint8,
                nodata=None,
            )
            save_raster(
                output_path=output_path,
                rasterio_profile=visualization_profile,
                raster_array=visualization_8bit_array,
                internal_nodata_mask_array=nodata_mask,
            )

        self._validate_save_horizon_visualization_parameters(
            dem_path=dem_path,
            overwrite=overwrite,
            save_svf_float_visualization=save_svf_float_visualization,
            output_svf_float_visualization_path=output_svf_float_visualization_path,
            save_svf_8bit_visualization=save_svf_8bit_visualization,
            output_svf_8bit_visualization_path=output_svf_8bit_visualization_path,
            save_asvf_float_visualization=save_asvf_float_visualization,
            output_asvf_float_visualization_path=output_asvf_float_visualization_path,
            save_asvf_8bit_visualization=save_asvf_8bit_visualization,
            output_asvf_8bit_visualization_path=output_asvf_8bit_visualization_path,
            save_opns_float_visualization=save_opns_float_visualization,
            output_opns_float_visualization_path=output_opns_float_visualization_path,
            save_opns_8bit_visualization=save_opns_8bit_visualization,
            output_opns_8bit_visualization_path=output_opns_8bit_visualization_path,
        )

        dem_file: rasterio.io.DatasetReader
        with rasterio.open(dem_path, "r") as dem_file:
            dem_transform = dem_file.profile["transform"]
            horizon_visualization_result = self.compute_horizon_visualizations(
                dem=dem_file.read(0),
                dem_resolution=(abs(dem_file.transform[0]) + abs(dem_file.transform[4]))
                / 2,
                dem_nodata=dem_file.nodata,
                vertical_exaggeration_factor=vertical_exaggeration_factor,
                compute_svf=save_svf_float_visualization or save_svf_8bit_visualization,
                compute_opns=save_opns_float_visualization
                or save_opns_8bit_visualization,
                compute_asvf=save_asvf_float_visualization
                or save_asvf_8bit_visualization,
            )
        svf_visualization_array = horizon_visualization_result.sky_view_factor
        asvf_visualization_array = (
            horizon_visualization_result.anisotropic_sky_view_factor
        )
        opns_visualization_array = horizon_visualization_result.openness
        if save_svf_float_visualization:
            _save_float_visualization(
                visualization_array=svf_visualization_array,
                output_path=output_svf_float_visualization_path,
            )
        if save_svf_8bit_visualization:
            assert svf_visualization_array is not None
            svf_visualization_8bit_array = (
                self.sky_view_factor.convert_visualization_to_8bit(
                    visualization_array=svf_visualization_array, no_data=np.nan
                )
            )
            _save_8bit_visualization(
                visualization_8bit_array=svf_visualization_8bit_array,
                output_path=output_svf_8bit_visualization_path,
            )
        if save_asvf_float_visualization:
            _save_float_visualization(
                visualization_array=asvf_visualization_array,
                output_path=output_asvf_float_visualization_path,
            )
        if save_asvf_8bit_visualization:
            assert asvf_visualization_array is not None
            asvf_visualization_8bit_array = (
                self.anisotropic_sky_view_factor.convert_visualization_to_8bit(
                    visualization_array=asvf_visualization_array, no_data=np.nan
                )
            )
            _save_8bit_visualization(
                visualization_8bit_array=asvf_visualization_8bit_array,
                output_path=output_asvf_8bit_visualization_path,
            )
        if save_opns_float_visualization:
            _save_float_visualization(
                visualization_array=opns_visualization_array,
                output_path=output_opns_float_visualization_path,
            )
        if save_opns_8bit_visualization:
            assert opns_visualization_array is not None
            opns_visualization_8bit_array = (
                self.positive_openness.convert_visualization_to_8bit(
                    visualization_array=opns_visualization_array, no_data=np.nan
                )
            )
            _save_8bit_visualization(
                visualization_8bit_array=opns_visualization_8bit_array,
                output_path=output_opns_8bit_visualization_path,
            )


class LocalDominance(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.LOCAL_DOMINANCE
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.50, maximum=1.80
    )
    _DEFAULT_MINIMUM_RADIUS: int = 10
    _DEFAULT_MAXIMUM_RADIUS: int = 20
    _DEFAULT_RADIAL_DISTANCE_STEP: int = 1
    _DEFAULT_ANGULAR_STEP: int = 15
    _DEFAULT_OBSERVER_HEIGHT: float = 1.7

    def __init__(
        self,
        minimum_radius: int = _DEFAULT_MINIMUM_RADIUS,
        maximum_radius: int = _DEFAULT_MAXIMUM_RADIUS,
        radial_distance_step: int = _DEFAULT_RADIAL_DISTANCE_STEP,
        angular_step: int = _DEFAULT_ANGULAR_STEP,
        observer_height: float = _DEFAULT_OBSERVER_HEIGHT,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.minimum_radius = minimum_radius
        self.maximum_radius = maximum_radius
        self.radial_distance_step = radial_distance_step
        self.angular_step = angular_step
        self.observer_height = observer_height
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_LD_R_M{self.minimum_radius}-{self.maximum_radius}_DI{self.radial_distance_step}_"
            f"A{self.angular_step}_OH{self.observer_height}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return local_dominance(
            dem=dem,
            minimum_radius=self.minimum_radius,
            maximum_radius=self.maximum_radius,
            radial_distance_step=self.radial_distance_step,
            angular_step=self.angular_step,
            observer_height=self.observer_height,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )


class SkyIllumination:  # (RVTVisualization)
    RVT_VISUALIZATION_NAME = RVTVisualizationName.SKY_ILLUMINATION
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.PERCENT, minimum=0.25, maximum=0.00
    )

    def __init__(self):
        raise NotImplemented  # TODO: Fix visualization and implement


class MultiScaleReliefModel(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.MULTI_SCALE_RELIEF_MODEL
    _DEFAULT_TO_8BIT: To8bit = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=-2.50, maximum=2.50
    )
    _DEFAULT_MINIMUM_FEATURE_SIZE: float = 0.0
    _DEFAULT_MAXIMUM_FEATURE_SIZE: float = 20.0
    _DEFAULT_SCALING_FACTOR: int = 2

    def __init__(
        self,
        minimum_feature_size: float = _DEFAULT_MINIMUM_FEATURE_SIZE,
        maximum_feature_size: float = _DEFAULT_MAXIMUM_FEATURE_SIZE,
        scaling_factor: int = _DEFAULT_SCALING_FACTOR,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.minimum_feature_size = minimum_feature_size
        self.maximum_feature_size = maximum_feature_size
        self.scaling_factor = scaling_factor
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_MSRM_FS{self.minimum_feature_size}-{self.maximum_feature_size}_"
            f"S{self.scaling_factor}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return multi_scale_relief_model(
            dem=dem,
            resolution=dem_resolution,
            minimum_feature_size=self.minimum_feature_size,
            maximum_feature_size=self.maximum_feature_size,
            scaling_factor=self.scaling_factor,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )


class MultiScaleTopographicPosition(RVTVisualization):
    RVT_VISUALIZATION_NAME = RVTVisualizationName.MULTI_SCALE_TOPOGRAPHIC_POSITION
    _DEFAULT_TO_8BIT = To8bit(
        normalization_mode=NormalizationMode.VALUE, minimum=0.00, maximum=1.00
    )
    _DEFAULT_LOCAL_SCALE: Tuple[int, int, int] = (3, 21, 2)
    _DEFAULT_MESO_SCALE: Tuple[int, int, int] = (23, 203, 18)
    _DEFAULT_BROAD_SCALE: Tuple[int, int, int] = (223, 2023, 180)
    _DEFAULT_LIGHTNESS: float = 1.2

    def __init__(
        self,
        local_scale: Tuple[int, int, int] = _DEFAULT_LOCAL_SCALE,
        meso_scale: Tuple[int, int, int] = _DEFAULT_MESO_SCALE,
        broad_scale: Tuple[int, int, int] = _DEFAULT_BROAD_SCALE,
        lightness: float = _DEFAULT_LIGHTNESS,
        to_8bit: To8bit = _DEFAULT_TO_8BIT,
    ):
        self.local_scale = local_scale
        self.meso_scale = meso_scale
        self.broad_scale = broad_scale
        self.lightness = lightness
        self.to_8bit = to_8bit

    def get_visualization_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        dem_file_name = dem_file_name.rstrip(".tif").rstrip(".tiff")
        output_name = (
            f"{dem_file_name}_MSTP_{self.local_scale[1]}_{self.meso_scale[1]}_{self.broad_scale[1]}_"
            f"L{self.lightness}"
        )
        if is_8bit:
            output_name += "_8bit"
        return output_name + ".tif"

    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        return multi_scale_topographic_position(
            dem=dem,
            local_scale=self.local_scale,
            meso_scale=self.meso_scale,
            broad_scale=self.broad_scale,
            lightness=self.lightness,
            vertical_exaggeration_factor=vertical_exaggeration_factor,
            no_data=dem_nodata,
        )


@dataclass
class RVTVisualizationFactory:
    """
    Class to store parameter values for all RVT visualizations, class also contain methods related to visualizations
    (methods like saving visualizations, getting visualizations file names, storing visualizations parameters to file,
     ...).
    """

    overwrite: bool = False
    """If True, overwrite visualizations when saving if they already exist."""
    vertical_exaggeration_factor: float = 1.0
    """Vertical exaggeration to apply to DEM before computing visualizations."""
    # Slope
    enable_slope: bool = False
    save_slope_float: bool = False
    save_slope_8bit: bool = False
    slope: Slope = Slope()
    # Shadow
    enable_shadow: bool = False
    save_shadow_float: bool = False
    save_shadow_8bit: bool = False
    shadow: Shadow = Shadow()
    # Hillshade
    enable_hillshade: bool = True
    save_hillshade_float: bool = True
    save_hillshade_8bit: bool = False
    hillshade: Hillshade = Hillshade()
    # Multiple Directions Hillshade
    enable_multiple_directions_hillshade: bool = False
    save_multiple_directions_hillshade_float: bool = False
    save_multiple_directions_hillshade_8bit: bool = False
    multiple_directions_hillshade: MultipleDirectionsHillshade = (
        MultipleDirectionsHillshade()
    )
    # Simple Local Relief Model
    enable_simple_local_relief_model: bool = False
    save_simple_local_relief_model_float: bool = False
    save_simple_local_relief_model_8bit: bool = False
    simple_local_relief_model: SimpleLocalReliefModel = SimpleLocalReliefModel()
    # Sky-View Factor
    enable_sky_view_factor: bool = False
    save_sky_view_factor_float: bool = False
    save_sky_view_factor_8bit: bool = False
    # Anisotropic Sky-View Factor
    enable_anisotropic_sky_view_factor: bool = False
    save_anisotropic_sky_view_factor_float: bool = False
    save_anisotropic_sky_view_factor_8bit: bool = False
    # Positive Openness
    enable_positive_openness: bool = False
    save_positive_openness_float: bool = False
    save_positive_openness_8bit: bool = False
    # Negative Openness
    enable_negative_openness: bool = False
    save_negative_openness_float: bool = False
    save_negative_openness_8bit: bool = False
    # Horizon Visualizations (Sky-View Factor, Anisotropic Sky-View factor, Openness)
    horizon_visualizations: HorizonVisualizations = HorizonVisualizations()
    # Sky Illumination # TODO: Fix
    # save_sky_illumination: bool = False
    # save_sky_illumination_float: bool = False
    # save_sky_illumination_8bit: bool = False
    # sky_illumination: SkyIllumination = SkyIllumination()
    # Local Dominance
    enable_local_dominance: bool = False
    save_local_dominance_float: bool = False
    save_local_dominance_8bit: bool = False
    local_dominance: LocalDominance = LocalDominance()
    # Multi Scale Relief Model
    enable_multi_scale_relief_model: bool = False
    save_multi_scale_relief_model_float: bool = False
    save_multi_scale_relief_model_8bit: bool = False
    multi_scale_relief_model: MultiScaleReliefModel = MultiScaleReliefModel()
    # Multi Scale Topographic Position
    enable_multi_scale_topographic_position: bool = False
    save_multi_scale_topographic_position_float: bool = False
    save_multi_scale_topographic_position_8bit: bool = False
    multi_scale_topographic_position: MultiScaleTopographicPosition = (
        MultiScaleTopographicPosition()
    )

    compute_multiple_tile_limit = 10000 * 10000
    """
    If number of pixels is bigger than defined value, visualizations are computed tile by tile.
    This is needed because some DEMs are too big (don't fit into memory) so output visualizations need to be computed
    part by part.
    """
    compute_tile_size = (
        4000,
        4000,
    )
    """Define the size of single tile if `compute_multiple_tile_limit` is exceeded."""

    def save_parameters_to_file(self, file_path: Path) -> None:
        """Save parameters to JSON file."""
        class_parameters_dict = deepcopy(self).__dict__
        for parameter, parameter_value in class_parameters_dict.items():
            if issubclass(type(parameter_value), RVTVisualization) or isinstance(
                parameter_value, HorizonVisualizations
            ):
                sub_class_parameters_dict = parameter_value.__dict__
                for (
                    sub_class_parameter,
                    sub_class_parameter_value,
                ) in sub_class_parameters_dict.items():
                    if isinstance(sub_class_parameter_value, RVTVisualizationName):
                        # remove RVTVisualizationName parameter since key already defines it
                        del sub_class_parameters_dict[sub_class_parameter]
                    elif isinstance(sub_class_parameter_value, To8bit):
                        sub_class_parameters_dict[
                            sub_class_parameter
                        ] = sub_class_parameter_value.to_string()
                    elif issubclass(type(sub_class_parameter_value), Enum):
                        sub_class_parameters_dict[
                            sub_class_parameter
                        ] = sub_class_parameter_value.value
                class_parameters_dict[parameter] = sub_class_parameters_dict

        with open(file_path, "w") as json_file:
            json.dump(class_parameters_dict, json_file, indent=4)

    @classmethod
    def read_parameters_from_file(cls, file_path: Path) -> RVTVisualizationFactory:
        """Reads parameters from JSON file."""

        def _init_visualization_class(
            visualization_parameters_dict: Dict[str, Any],
            visualization_type: Type[RVTVisualization],
        ) -> RVTVisualization:
            for (
                visualization_attribute_name,
                visualization_attribute_value,
            ) in visualization_parameters_dict.items():
                visualization_attribute_type = type(
                    # Initialize visualization_type to check the type of default values.
                    getattr(visualization_type(), visualization_attribute_name)
                )
                if visualization_attribute_name == "to_8bit":
                    visualization_parameters_dict["to_8bit"] = To8bit.from_string(
                        string=visualization_parameters_dict["to_8bit"]
                    )
                else:
                    # Default values in child classes of `RVTVisualization` should match their annotation
                    # (i.e. _DEFAULT_SUN_ELEVATION should be float number not integer (35.0 not 35) in order to get the
                    # right type when reading parameters from file. Currently, attributes of type `RVTVisualization` in
                    # `RVTVisualizationFactory` are first initialized with default values to get the right type for
                    # every parameter. This type is then used to convert string value to the right type.
                    # This issue could be solved with PyDantic BaseModel and Field, for now we didn't want to introduce
                    # dependency to PyDantic.
                    visualization_parameters_dict[
                        visualization_attribute_name
                    ] = visualization_attribute_type(
                        visualization_parameters_dict[visualization_attribute_name]
                    )

            return visualization_type(**visualization_parameters_dict)

        with open(file_path, "r") as json_file:
            parameters_dict = json.load(json_file)

        for attribute_name in cls.__annotations__:
            attribute_type = type(getattr(cls, attribute_name))
            if issubclass(attribute_type, RVTVisualization) or issubclass(
                attribute_type, HorizonVisualizations
            ):
                if attribute_name not in parameters_dict:
                    continue
                parameters_dict[attribute_name] = _init_visualization_class(
                    visualization_parameters_dict=parameters_dict[attribute_name],
                    visualization_type=attribute_type,
                )
        return cls(**parameters_dict)

    def _get_visualization(
        self,
        visualization: RVTVisualization,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        visualization_array = visualization.compute_visualization(
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            vertical_exaggeration_factor=self.vertical_exaggeration_factor,
        )
        if is_8bit:
            return visualization.convert_visualization_to_8bit(
                visualization_array=visualization_array
            )
        return visualization_array

    def get_shadow_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return self.shadow.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_shadow_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.shadow.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_shadow_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.shadow,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_slope_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return self.slope.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_slope_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.slope.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_slope_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.slope,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_hillshade_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return self.hillshade.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_hillshade_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.hillshade.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_hillshade_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.hillshade,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_multiple_directions_hillshade_file_name(
        self, dem_file_name: str, is_8bit: bool
    ) -> str:
        return self.multiple_directions_hillshade.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_multiple_directions_hillshade_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.multiple_directions_hillshade.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_multiple_directions_hillshade_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        if is_8bit:
            return self.multiple_directions_hillshade.compute_8bit_visualization(
                dem=dem,
                dem_resolution=dem_resolution,
                dem_nodata=dem_nodata,
                vertical_exaggeration_factor=self.vertical_exaggeration_factor,
            )
        return self.multiple_directions_hillshade.compute_visualization(
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            vertical_exaggeration_factor=self.vertical_exaggeration_factor,
        )

    def get_simple_local_relief_model_file_name(
        self, dem_file_name: str, is_8bit: bool
    ) -> str:
        return self.simple_local_relief_model.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_simple_local_relief_model_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.simple_local_relief_model.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_simple_local_relief_model_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.simple_local_relief_model,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_sky_view_factor_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return self.horizon_visualizations.sky_view_factor.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_sky_view_factor_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.horizon_visualizations.sky_view_factor.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_sky_view_factor_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.horizon_visualizations.sky_view_factor,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_anisotropic_sky_view_factor_file_name(
        self, dem_file_name: str, is_8bit: bool
    ) -> str:
        return self.horizon_visualizations.anisotropic_sky_view_factor.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_anisotropic_sky_view_factor_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.horizon_visualizations.anisotropic_sky_view_factor.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_anisotropic_sky_view_factor_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.horizon_visualizations.anisotropic_sky_view_factor,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_positive_openness_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return (
            self.horizon_visualizations.positive_openness.get_visualization_file_name(
                dem_file_name=dem_file_name, is_8bit=is_8bit
            )
        )

    def get_positive_openness_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.horizon_visualizations.positive_openness.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_positive_openness_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.horizon_visualizations.positive_openness,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_negative_openness_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return (
            self.horizon_visualizations.negative_openness.get_visualization_file_name(
                dem_file_name=dem_file_name, is_8bit=is_8bit
            )
        )

    def get_negative_openness_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.horizon_visualizations.negative_openness.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_negative_openness_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.horizon_visualizations.negative_openness,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_sky_illumination_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        raise NotImplemented

    def get_sky_illumination_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        raise NotImplemented

    def get_sky_illumination_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        raise NotImplemented

    def get_local_dominance_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        return self.local_dominance.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_local_dominance_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.local_dominance.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_local_dominance_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.local_dominance,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_multi_scale_relief_model_file_name(
        self, dem_file_name: str, is_8bit: bool
    ) -> str:
        return self.multi_scale_relief_model.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_multi_scale_relief_model_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.multi_scale_relief_model.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_multi_scale_relief_model_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.multi_scale_relief_model,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_multi_scale_topographic_position_file_name(
        self, dem_file_name: str, is_8bit: bool
    ) -> str:
        return self.multi_scale_topographic_position.get_visualization_file_name(
            dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_multi_scale_topographic_position_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
        return self.multi_scale_topographic_position.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )

    def get_multi_scale_topographic_position_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        is_8bit: bool = False,
    ) -> npt.NDArray[Any]:
        self._get_visualization(
            visualization=self.multi_scale_topographic_position,
            dem=dem,
            dem_resolution=dem_resolution,
            dem_nodata=dem_nodata,
            is_8bit=is_8bit,
        )

    def get_visualization_file_name(
        self,
        visualization_name: RVTVisualizationName,
        dem_file_name: str,
        is_8bit: bool,
    ) -> str:
        for field in fields(self.__class__):
            field_value = getattr(self, field.name)
            if issubclass(field_value, RVTVisualization):
                if field_value.RVT_VISUALIZATION_NAME == visualization_name:
                    return field_value.get_visualization_file_name(
                        dem_file_name=dem_file_name, is_8bit=is_8bit
                    )
        raise ValueError(f"Unknown visualization name {visualization_name}!")

    def get_visualization_path(
        self,
        visualization_name: RVTVisualizationName,
        directory_path: Path,
        dem_file_name: str,
        is_8bit: bool,
    ) -> Path:
        # Loop through the fields of RVTVisualizationFactory.
        for field in fields(self.__class__):
            field_value = getattr(self, field.name)
            # If field value is child class of RVTVisualization.
            if issubclass(type(field_value), RVTVisualization):
                if field_value.RVT_VISUALIZATION_NAME == visualization_name:
                    return field_value.get_visualization_path(
                        directory_path=directory_path,
                        dem_file_name=dem_file_name,
                        is_8bit=is_8bit,
                    )
            # If field value is HorizonVisualizations.
            if isinstance(field_value, HorizonVisualizations):
                # Loop through the fields of HorizonVisualizations.
                for horizon_vis_attr_name in dir(field_value):
                    horizon_vis_field_value = getattr(
                        field_value, horizon_vis_attr_name
                    )
                    # If horizon visualization field value is child class of RVTVisualization.
                    if issubclass(type(horizon_vis_field_value), RVTVisualization):
                        if (
                            horizon_vis_field_value.RVT_VISUALIZATION_NAME
                            == visualization_name
                        ):
                            return horizon_vis_field_value.get_visualization_path(
                                directory_path=directory_path,
                                dem_file_name=dem_file_name,
                                is_8bit=is_8bit,
                            )
        raise ValueError(f"Unknown visualization name {visualization_name}!")

    def save_visualizations(self) -> None:
        pass
        # TODO: ZiM you stayed here


#  TODO: Rewrite and clean methods/functions bellow
#
#     def create_log_file(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         compute_time: Optional[float] = None,
#     ) -> None:
#         """
#         Creates log file in custom_dir, if custom_dir=None it creates it in dem directory (dem_path).
#         Be aware, all default parameters have to be right! Parameter compute_time is in seconds.
#         """
#         dict_arr_res = get_raster_arr(raster_path=dem_path)
#         resolution = dict_arr_res["resolution"]
#         arr_shape = np.array(dict_arr_res["array"]).shape
#         del dict_arr_res
#         nr_bands = 0
#         nr_cols = 0
#         nr_rows = 0
#         if len(arr_shape) == 3:
#             nr_bands = arr_shape[0]
#             nr_rows = arr_shape[1]
#             nr_cols = arr_shape[2]
#         elif len(arr_shape) == 2:
#             nr_bands = 1
#             nr_rows = arr_shape[0]
#             nr_cols = arr_shape[1]
#         dem_dir = Path(os.path.dirname(dem_path))
#         log_dir = dem_dir
#         if custom_dir is not None:
#             log_dir = custom_dir
#         dem_name = os.path.splitext(os.path.basename(dem_path))[0]
#         log_file_time = datetime.datetime.now()
#         log_file_time_str = log_file_time.strftime("%Y-%m-%d_%H-%M-%S")
#         log_name = "{}_vis_log_{}.txt".format(dem_name, log_file_time_str)
#         log_path = os.path.join(log_dir, log_name)
#         dat = open(log_path, "w")
#         dat.write(
#             "===============================================================================================\n"
#             "Relief RVTVisualization Toolbox (python), visualizations log\n"
#             "Copyright:\n"
#             "\tResearch Centre of the Slovenian Academy of Sciences and Arts\n"
#             "\tUniversity of Ljubljana, Faculty of Civil and Geodetic Engineering\n"
#             "===============================================================================================\n"
#         )
#         dat.write("\n\n\n")
#
#         dat.write(
#             "Processing info about visualizations\n"
#             "===============================================================================================\n\n"
#         )
#         dat.write("# Metadata of the input file\n\n")
#         dat.write("\tInput filename:\t\t{}\n".format(dem_path))
#         dat.write("\tNumber of rows:\t\t{}\n".format(nr_rows))
#         dat.write("\tNumber of columns:\t{}\n".format(nr_cols))
#         dat.write("\tNumber of bands:\t{}\n".format(nr_bands))
#         dat.write("\tResolution (x, y):\t{}, {}\n".format(resolution[0], resolution[1]))
#         dat.write("\n")
#
#         dat.write("# Selected visualization parameters\n")
#         dat.write("\tOverwrite: {}\n".format(self.overwrite))
#         dat.write(
#             "\tVertical exaggeration factor: {}\n".format(
#                 self.vertical_exaggeration_factor
#             )
#         )
#         if nr_rows * nr_cols > self.compute_multiple_tile_limit:
#             dat.write("\tCalculating tile by tile: {}\n".format("ON"))
#             dat.write(
#                 "\t\tTile block size: {}x{}\n".format(
#                     self.tile_size[0], self.tile_size[1]
#                 )
#             )
#         else:
#             dat.write("\tCalculating tile by tile: {}\n".format("OFF"))
#
#         dat.write("\n")
#
#         dat.write("# The following visualizations have been preformed:\n\n")
#         if self.hs_compute:
#             dat.write("\tHillshade\n")
#             dat.write("\t\ths_sun_el=\t\t{}\n".format(self.hs_sun_el))
#             dat.write("\t\ths_sun_azi=\t\t{}\n".format(self.hs_sun_azi))
#             if self.hs_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_hillshade_file_name(dem_path)
#                             )
#                         )
#                     )
#                 )
#             if self.hs_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\ths_bytscl=\t\t({}, {}, {})\n".format(
#                         self.hs_bytscl[0], self.hs_bytscl[1], self.hs_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir,
#                                 self.get_hillshade_file_name(dem_path, is_8bit=True),
#                             )
#                         )
#                     )
#                 )
#             if self.hs_shadow:
#                 dat.write("\t\t>> Output shadow file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_shadow_file_name(dem_path))
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.mhs_compute:
#             dat.write("\tMultiple directions hillshade\n")
#             dat.write("\t\tmhs_sun_el=\t\t{}\n".format(self.mhs_sun_el))
#             dat.write("\t\tmhs_nr_dir=\t\t{}\n".format(self.mhs_nr_dir))
#             if self.mhs_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_multi_hillshade_file_name(dem_path)
#                             )
#                         )
#                     )
#                 )
#             if self.mhs_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tmhs_bytscl=\t\t({}, {}, {})\n".format(
#                         self.mhs_bytscl[0], self.mhs_bytscl[1], self.mhs_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir,
#                                 self.get_multi_hillshade_file_name(
#                                     dem_path, is_8bit=True
#                                 ),
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.slp_compute:
#             dat.write("\tSlope gradient\n")
#             dat.write("\t\tslp_output_units=\t\t{}\n".format(self.slp_output_units))
#             if self.slp_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_slope_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.slp_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tslp_bytscl=\t\t({}, {}, {})\n".format(
#                         self.slp_bytscl[0], self.slp_bytscl[1], self.slp_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir,
#                                 self.get_slope_file_name(dem_path, is_8bit=True),
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.slrm_compute:
#             dat.write("\tSimple local relief model\n")
#             dat.write("\t\tslrm_rad_cell=\t\t{}\n".format(self.slrm_rad_cell))
#             if self.slrm_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_slrm_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.slrm_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tslrm_bytscl=\t\t({}, {}, {})\n".format(
#                         self.slrm_bytscl[0], self.slrm_bytscl[1], self.slrm_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_slrm_file_name(dem_path, is_8bit=True)
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.msrm_compute:
#             dat.write("\tMulti-scale relief model\n")
#             dat.write("\t\tmsrm_feature_min=\t\t{}\n".format(self.msrm_feature_min))
#             dat.write("\t\tmsrm_feature_max=\t\t{}\n".format(self.msrm_feature_max))
#             dat.write(
#                 "\t\tmsrm_scaling_factor=\t\t{}\n".format(self.msrm_scaling_factor)
#             )
#             if self.msrm_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_msrm_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.msrm_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tmsrm_bytscl=\t\t({}, {}, {})\n".format(
#                         self.msrm_bytscl[0], self.msrm_bytscl[1], self.msrm_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_msrm_file_name(dem_path, is_8bit=True)
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.svf_compute:
#             dat.write("\tSky-View Factor\n")
#             dat.write("\t\tnumber_of_directions=\t\t{}\n".format(self.svf_n_dir))
#             dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
#             dat.write("\t\tmaximum_search_radius=\t\t{}\n".format(self.svf_r_max))
#             if self.svf_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_svf_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.svf_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tsvf_bytscl=\t\t({}, {}, {})\n".format(
#                         self.svf_bytscl[0], self.svf_bytscl[1], self.svf_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_svf_file_name(dem_path, is_8bit=True)
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.asvf_compute:
#             dat.write("\tAnisotropic Sky-View Factor\n")
#             dat.write("\t\tnumber_of_directions=\t\t{}\n".format(self.svf_n_dir))
#             dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
#             dat.write("\t\tmaximum_search_radius=\t\t{}\n".format(self.svf_r_max))
#             dat.write("\t\tanisotropy_level=\t\t{}\n".format(self.asvf_level))
#             dat.write("\t\tdirection_of_anisotropy=\t\t{}\n".format(self.asvf_dir))
#             if self.svf_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_asvf_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.svf_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tasvf_bytscl=\t\t({}, {}, {})\n".format(
#                         self.asvf_bytscl[0], self.asvf_bytscl[1], self.asvf_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_asvf_file_name(dem_path, is_8bit=True)
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.pos_opns_compute:
#             dat.write("\tOpenness - Positive\n")
#             dat.write("\t\tnumber_of_directions=\t\t{}\n".format(self.svf_n_dir))
#             dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
#             dat.write("\t\tmaximum_search_radius=\t\t{}\n".format(self.svf_r_max))
#             if self.svf_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_opns_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.svf_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tpos_opns_bytscl=\t\t({}, {}, {})\n".format(
#                         self.pos_opns_bytscl[0],
#                         self.pos_opns_bytscl[1],
#                         self.pos_opns_bytscl[2],
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_opns_file_name(dem_path, is_8bit=True)
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.neg_opns_compute:
#             dat.write("\tOpenness - Negative\n")
#             dat.write("\t\tnumber_of_directions=\t\t{}\n".format(self.svf_n_dir))
#             dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
#             dat.write("\t\tmaximum_search_radius=\t\t{}\n".format(self.svf_r_max))
#             if self.neg_opns_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(log_dir, self.get_neg_opns_file_name(dem_path))
#                         )
#                     )
#                 )
#             if self.neg_opns_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tneg_opns_bytscl=\t\t({}, {}, {})\n".format(
#                         self.neg_opns_bytscl[0],
#                         self.neg_opns_bytscl[1],
#                         self.neg_opns_bytscl[2],
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir,
#                                 self.get_neg_opns_file_name(dem_path, is_8bit=True),
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.sim_compute:
#             dat.write("\tSky illumination\n")
#             dat.write("\t\tsim_sky_mod=\t\t{}\n".format(self.sim_sky_mod))
#             dat.write("\t\tsim_compute_shadow=\t\t{}\n".format(self.sim_compute_shadow))
#             dat.write("\t\tsim_shadow_az=\t\t{}\n".format(self.sim_shadow_az))
#             dat.write("\t\tsim_shadow_el=\t\t{}\n".format(self.sim_shadow_el))
#             dat.write("\t\tsim_nr_dir=\t\t{}\n".format(self.sim_nr_dir))
#             dat.write("\t\tsim_shadow_dist=\t\t{}\n".format(self.sim_shadow_dist))
#             if self.sim_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_sky_illumination_file_name(dem_path)
#                             )
#                         )
#                     )
#                 )
#             if self.sim_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tsim_bytscl=\t\t({}, {}, {})\n".format(
#                         self.sim_bytscl[0], self.sim_bytscl[1], self.sim_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir,
#                                 self.get_sky_illumination_file_name(
#                                     dem_path, is_8bit=True
#                                 ),
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.ld_compute:
#             dat.write("\tLocal dominance\n")
#             dat.write("\t\tld_rad_inc=\t\t{}\n".format(self.ld_rad_inc))
#             dat.write("\t\tld_min_rad=\t\t{}\n".format(self.ld_min_rad))
#             dat.write("\t\tld_max_rad=\t\t{}\n".format(self.ld_max_rad))
#             dat.write("\t\tld_anglr_res=\t\t{}\n".format(self.ld_anglr_res))
#             dat.write("\t\tld_observer_h=\t\t{}\n".format(self.ld_observer_h))
#             if self.ld_save_float:
#                 dat.write("\t\t>> Output file:\n")
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir, self.get_local_dominance_file_name(dem_path)
#                             )
#                         )
#                     )
#                 )
#             if self.ld_save_8bit:
#                 dat.write("\t\t>> Output 8bit file:\n")
#                 dat.write(
#                     "\t\tld_bytscl=\t\t({}, {}, {})\n".format(
#                         self.ld_bytscl[0], self.ld_bytscl[1], self.ld_bytscl[2]
#                     )
#                 )
#                 dat.write(
#                     "\t\t\t{}\n".format(
#                         os.path.abspath(
#                             os.path.join(
#                                 log_dir,
#                                 self.get_local_dominance_file_name(
#                                     dem_path, is_8bit=True
#                                 ),
#                             )
#                         )
#                     )
#                 )
#             dat.write("\n")
#         if self.mstp_compute:
#             dat.write("\tMulti-scale topographic position\n")
#             dat.write(
#                 "\t\tmstp_local_scale=\t({}, {}, {})\n".format(
#                     self.mstp_local_scale[0],
#                     self.mstp_local_scale[1],
#                     self.mstp_local_scale[2],
#                 )
#             )
#             dat.write(
#                 "\t\tmstp_meso_scale=\t({}, {}, {})\n".format(
#                     self.mstp_meso_scale[0],
#                     self.mstp_meso_scale[1],
#                     self.mstp_meso_scale[2],
#                 )
#             )
#             dat.write(
#                 "\t\tmstp_broad_scale=\t({}, {}, {})\n".format(
#                     self.mstp_broad_scale[0],
#                     self.mstp_broad_scale[1],
#                     self.mstp_broad_scale[2],
#                 )
#             )
#             dat.write("\t\tmstp_lightness=\t\t{}\n".format(self.mstp_lightness))
#             dat.write("\t\t>> Output file:\n")
#             dat.write(
#                 "\t\t\t{}\n".format(
#                     os.path.abspath(
#                         os.path.join(log_dir, self.get_mstp_file_name(dem_path))
#                     )
#                 )
#             )
#             dat.write("\n")
#
#         if compute_time is not None:
#             dat.write("# Computation time: {:.3f}s".format(compute_time))
#         dat.close()
#
#
# def get_raster_arr(raster_path: Path) -> Dict[str, Any]:
#     # TODO: Refactor and replace with rasterio
#     """
#     Reads raster from raster_path and returns its array(value) and resolution.
#
#     Parameters
#     ----------
#     raster_path : str
#         Path to raster
#
#     Returns
#     -------
#     dict_out : dict
#         Returns {"array": array, "resolution": (x_res, y_res), "no_data": no_data} : dict("array": np.array,
#         "resolution": tuple(float, float), "no_data": float).
#         Returns dictionary with keys: array, resolution and no_data. Key resolution is tuple where first element is x
#         resolution and second is y resolution. Key no_data represent value of no data.
#     """
#     data_set = gdal.Open(raster_path)
#     gt = data_set.GetGeoTransform()
#     x_res = abs(gt[1])
#     y_res = abs(-gt[5])
#     bands = []
#     no_data = data_set.GetRasterBand(
#         1
#     ).GetNoDataValue()  # we assume that all the bands have same no_data val
#     if data_set.RasterCount == 1:  # only one band_number
#         array = np.array(data_set.GetRasterBand(1).ReadAsArray())
#         data_set = None
#         return {"array": array, "resolution": (x_res, y_res), "no_data": no_data}
#     else:  # multiple bands
#         for i_band in range(data_set.RasterCount):
#             i_band += 1
#             band = np.array(data_set.GetRasterBand(i_band).ReadAsArray())
#             if band is None:
#                 continue
#             else:
#                 bands.append(band)
#         data_set = None  # close dataset
#         return {
#             "array": np.array(bands),
#             "resolution": (x_res, y_res),
#             "no_data": no_data,
#         }
#
