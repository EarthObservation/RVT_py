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

from rvt.blender_functions import Normalization
from rvt.enums import (
    RVTVisualizationName,
    OpennessType,
    SlopeOutputUnit,
    SvfNoiseRemove,
    AnisotropyLevel,
    NormalizationMode,
)
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


class RVTVisualization(ABC):
    RVT_VISUALIZATION_NAME: RVTVisualizationName
    to_8bit: To8bit

    @abstractmethod
    def compute_visualization(
        self,
        dem: npt.NDArray[Any],
        dem_resolution: float,
        dem_nodata: Optional[float],
        vertical_exaggeration_factor: float,
    ) -> npt.NDArray[Any]:
        pass

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

    def compute_visualization_from_dem_file(self, dem_path: Path) -> npt.NDArray[Any]:
        pass  # TODO: implement

    def save_visualization(
        self,
        dem_path: Path,
        output_visualization_path: Path,
        output_visualization: npt.NDArray,
    ) -> None:
        pass  # TODO: implement


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
            f"{dem_file_name}_shadow_A{self.sun_azimuth}_H{self.sun_elevation}"
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
            f"{dem_file_name}_SVF_R{self.maximum_search_radius}_D{self.number_of_directions}_"
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
        if OpennessType.POSITIVE:
            return RVTVisualizationName.POSITIVE_OPENNESS
        elif OpennessType.NEGATIVE:
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
            f"{dem_file_name}_MSRM_F_M{self.minimum_feature_size}-{self.maximum_feature_size}_"
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
    local_dominance = LocalDominance()
    # Multi Scale Relief Model
    enable_multi_scale_relief_model: bool = False
    save_multi_scale_relief_model_float: bool = False
    save_multi_scale_relief_model_8bit: bool = False
    multi_scale_relief_model = MultiScaleReliefModel()
    # Multi Scale Topographic Position
    enable_multi_scale_topographic_position: bool = False
    save_multi_scale_topographic_position_float: bool = False
    save_multi_scale_topographic_position_8bit: bool = False
    multi_scale_topographic_position = MultiScaleTopographicPosition()

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

    def get_sky_illumination_file_name(self, dem_file_name: str, is_8bit: bool) -> str:
        raise NotImplemented

    def get_sky_illumination_path(
        self, directory_path: Path, dem_file_name: str, is_8bit: bool
    ) -> Path:
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
        for field in fields(self.__class__):
            field_value = getattr(self, field.name)
            if issubclass(field_value, RVTVisualization):
                if field_value.RVT_VISUALIZATION_NAME == visualization_name:
                    return field_value.get_visualization_path(
                        directory_path=directory_path,
                        dem_file_name=dem_file_name,
                        is_8bit=is_8bit,
                    )
        raise ValueError(f"Unknown visualization name {visualization_name}!")


#  TODO: Rewrite and clean methods/functions bellow

#     def float_to_8bit(
#         self,
#         float_arr: npt.NDArray[Any],
#         visualization: RVTVisualizationName,
#         x_res: Optional[float] = None,
#         y_res: Optional[float] = None,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         """
#         Converts (byte scale) float visualization to 8bit. Resolution (x_res, y_res) and no_data needed only for
#         multiple directions hillshade! Method first normalize then byte scale (0-255).
#         """
#         if visualization == RVTVisualizationName.HILLSHADE:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="hs",
#                 image=float_arr,
#                 min_norm=self.hs_bytscl[1],
#                 max_norm=self.hs_bytscl[2],
#                 normalization=self.hs_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.SLOPE:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="slp",
#                 image=float_arr,
#                 min_norm=self.slp_bytscl[1],
#                 max_norm=self.slp_bytscl[2],
#                 normalization=self.slp_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.SHADOW:
#             return float_arr
#         elif visualization == RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE:
#             if x_res is None or y_res is None:
#                 raise ValueError("Resolution (x_res, y_res) needs to be specified!")
#             # Be careful when multihillshade we input dem, because we have to calculate hillshade in 3 directions
#             red_band_arr = rvt.vis.hillshade(
#                 dem=float_arr,
#                 resolution_x=x_res,
#                 resolution_y=y_res,
#                 sun_elevation=self.mhs_sun_el,
#                 sun_azimuth=315,
#                 no_data=no_data,
#             )
#             green_band_arr = rvt.vis.hillshade(
#                 dem=float_arr,
#                 resolution_x=x_res,
#                 resolution_y=y_res,
#                 sun_elevation=self.mhs_sun_el,
#                 sun_azimuth=22.5,
#                 no_data=no_data,
#             )
#             blue_band_arr = rvt.vis.hillshade(
#                 dem=float_arr,
#                 resolution_x=x_res,
#                 resolution_y=y_res,
#                 sun_elevation=self.mhs_sun_el,
#                 sun_azimuth=90,
#                 no_data=no_data,
#             )
#             if (
#                 self.mhs_bytscl[0].lower() == "percent"
#                 or self.slp_bytscl[0].lower() == "perc"
#             ):
#                 red_band_arr = rvt.blend_func.normalize_perc(
#                     image=red_band_arr,
#                     minimum=self.mhs_bytscl[1],
#                     maximum=self.mhs_bytscl[2],
#                 )
#                 red_band_arr = rvt.vis.byte_scale(
#                     data=red_band_arr, no_data=np.nan, c_min=0, c_max=1
#                 )
#                 green_band_arr = rvt.blend_func.normalize_perc(
#                     image=green_band_arr,
#                     minimum=self.mhs_bytscl[1],
#                     maximum=self.mhs_bytscl[2],
#                 )
#                 green_band_arr = rvt.vis.byte_scale(
#                     data=green_band_arr, no_data=np.nan, c_min=0, c_max=1
#                 )
#                 blue_band_arr = rvt.blend_func.normalize_perc(
#                     image=blue_band_arr,
#                     minimum=self.mhs_bytscl[1],
#                     maximum=self.mhs_bytscl[2],
#                 )
#                 blue_band_arr = rvt.vis.byte_scale(
#                     data=blue_band_arr, no_data=np.nan, c_min=0, c_max=1
#                 )
#             else:  # self.mhs_bytscl[0] == "value"
#                 red_band_arr = rvt.blend_func.normalize_lin(
#                     image=red_band_arr,
#                     minimum=self.mhs_bytscl[1],
#                     maximum=self.mhs_bytscl[2],
#                 )
#                 red_band_arr = rvt.vis.byte_scale(
#                     data=red_band_arr, no_data=np.nan, c_min=0, c_max=1
#                 )
#                 green_band_arr = rvt.blend_func.normalize_lin(
#                     image=green_band_arr,
#                     minimum=self.mhs_bytscl[1],
#                     maximum=self.mhs_bytscl[2],
#                 )
#                 green_band_arr = rvt.vis.byte_scale(
#                     data=green_band_arr, no_data=np.nan, c_min=0, c_max=1
#                 )
#                 blue_band_arr = rvt.blend_func.normalize_lin(
#                     image=blue_band_arr,
#                     minimum=self.mhs_bytscl[1],
#                     maximum=self.mhs_bytscl[2],
#                 )
#                 blue_band_arr = rvt.vis.byte_scale(
#                     data=blue_band_arr, no_data=np.nan, c_min=0, c_max=1
#                 )
#             multi_hillshade_8bit_arr = np.array(
#                 [red_band_arr, green_band_arr, blue_band_arr]
#             )
#             return multi_hillshade_8bit_arr
#         elif visualization == RVTVisualizationName.SIMPLE_LOCAL_RELIEF_MODEL:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="simple_local_relief_model",
#                 image=float_arr,
#                 min_norm=self.slrm_bytscl[1],
#                 max_norm=self.slrm_bytscl[2],
#                 normalization=self.slrm_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.SKY_VIEW_FACTOR:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="svf",
#                 image=float_arr,
#                 min_norm=self.svf_bytscl[1],
#                 max_norm=self.svf_bytscl[2],
#                 normalization=self.svf_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.ANISOTROPIC_SKY_VIEW_FACTOR:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="asvf",
#                 image=float_arr,
#                 min_norm=self.asvf_bytscl[1],
#                 max_norm=self.asvf_bytscl[2],
#                 normalization=self.asvf_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.POSITIVE_OPENNESS:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="pos_opns",
#                 image=float_arr,
#                 min_norm=self.pos_opns_bytscl[1],
#                 max_norm=self.pos_opns_bytscl[2],
#                 normalization=self.pos_opns_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.NEGATIVE_OPENNESS:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="neg_opns",
#                 image=float_arr,
#                 min_norm=self.neg_opns_bytscl[1],
#                 max_norm=self.neg_opns_bytscl[2],
#                 normalization=self.neg_opns_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.SKY_ILLUMINATION:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="sim",
#                 image=float_arr,
#                 min_norm=self.sim_bytscl[1],
#                 max_norm=self.sim_bytscl[2],
#                 normalization=self.sim_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.LOCAL_DOMINANCE:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="ld",
#                 image=float_arr,
#                 min_norm=self.ld_bytscl[1],
#                 max_norm=self.ld_bytscl[2],
#                 normalization=self.ld_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.MULTI_SCALE_RELIEF_MODEL:
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="multi_scale_relief_model",
#                 image=float_arr,
#                 min_norm=self.msrm_bytscl[1],
#                 max_norm=self.msrm_bytscl[2],
#                 normalization=self.msrm_bytscl[0],
#             )
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         elif visualization == RVTVisualizationName.MULTI_SCALE_TOPOGRAPHIC_POSITION:
#             # This might not be necessary, as all multi_scale_topographic_position data should already be between 0 and 1
#             norm_arr = rvt.blend_func.normalize_image(
#                 visualization="multi_scale_topographic_position",
#                 image=float_arr,
#                 min_norm=self.mstp_bytscl[1],
#                 max_norm=self.mstp_bytscl[2],
#                 normalization=self.mstp_bytscl[0],
#             )
#
#             return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
#         else:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.float_to_8bit: Wrong visualization (visualization) parameter!"
#             )
#
#     def get_slope(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution_x: float,
#         resolution_y: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         slope_arr = rvt.vis.slope_aspect(
#             dem=dem_arr,
#             resolution_x=resolution_x,
#             resolution_y=resolution_y,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             output_unit=self.slp_output_units,
#             no_data=no_data,
#         )["slope"]
#         return slope_arr
#
#     def save_slope(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """Calculates and saves Slope from dem (dem_path) with default parameters. If custom_dir is None it saves
#         in dem directory else in custom_dir. If path to file already exists we can overwrite file (overwrite=0) or
#         not (overwrite=1). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255)."""
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.slp_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.slp_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_slope: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_slope: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             slope_path = self.get_slope_path(dem_path)
#             slope_8bit_path = self.get_slope_path(dem_path, is_8bit=True)
#         else:
#             slope_path = custom_dir / self.get_slope_file_name(dem_path)
#             slope_8bit_path = custom_dir / self.get_slope_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(slope_8bit_path)
#                 and os.path.isfile(slope_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(slope_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(slope_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if (
#             dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit
#         ):  # tile by tile calculation
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.SLOPE,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#             slope_arr = self.get_slope(
#                 dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res, no_data=no_data
#             )
#             if save_float:
#                 if (
#                     os.path.isfile(slope_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=slope_path,
#                         out_raster_arr=slope_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(slope_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     slope_8bit_arr = self.float_to_8bit(
#                         float_arr=slope_arr, visualization=RVTVisualizationName.SLOPE
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=slope_8bit_path,
#                         out_raster_arr=slope_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_shadow(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         shadow_arr = rvt.vis.shadow_horizon(
#             dem=dem_arr,
#             resolution=resolution,
#             shadow_az=self.hs_sun_azi,
#             shadow_el=self.hs_sun_el,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )["shadow"]
#         return shadow_arr
#
#     def get_hillshade(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution_x: float,
#         resolution_y: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         hillshade_arr = rvt.vis.hillshade(
#             dem=dem_arr,
#             resolution_x=resolution_x,
#             resolution_y=resolution_y,
#             sun_azimuth=self.hs_sun_azi,
#             sun_elevation=self.hs_sun_el,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return hillshade_arr
#
#     def save_hillshade(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#         save_shadow: Optional[bool] = None,
#     ) -> bool:
#         """Calculates and saves Hillshade from dem (dem_path) with default parameters. If custom_dir is None it saves
#         in dem directory else in custom_dir. If path to file already exists we can overwrite file (overwrite=1)
#         or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255)."""
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.hs_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.hs_save_8bit)
#         # if save_shadow is None it takes boolean from default (self)
#         if save_shadow is None:
#             save_shadow = bool(self.hs_shadow)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_hillshade: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_hillshade: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             hillshade_path = self.get_hillshade_path(dem_path)
#             hillshade_8bit_path = self.get_hillshade_path(dem_path, is_8bit=True)
#             shadow_path = self.get_shadow_path(dem_path)
#         else:
#             hillshade_path = custom_dir / self.get_hillshade_file_name(dem_path)
#             hillshade_8bit_path = custom_dir / self.get_hillshade_file_name(
#                 dem_path, is_8bit=True
#             )
#             shadow_path = custom_dir / self.get_shadow_file_name(dem_path)
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit and save_shadow:
#             if (
#                 os.path.isfile(hillshade_8bit_path)
#                 and os.path.isfile(hillshade_path)
#                 and os.path.isfile(shadow_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit and not save_shadow:
#             if os.path.isfile(hillshade_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit and not save_shadow:
#             if os.path.isfile(hillshade_8bit_path) and not self.overwrite:
#                 return False
#         elif save_float and not save_8bit and save_shadow:
#             if (
#                 os.path.isfile(hillshade_path)
#                 and os.path.isfile(shadow_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif not save_float and not save_8bit and save_shadow:
#             if (
#                 os.path.isfile(hillshade_8bit_path)
#                 and os.path.isfile(shadow_path)
#                 and not self.overwrite
#             ):
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.HILLSHADE,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             if save_shadow:
#                 rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                     RVT_VISUALIZATION_NAME=RVTVisualizationName.SHADOW,
#                     rvt_default=self,
#                     dem_path=Path(dem_path),
#                     output_dir_path=Path(custom_dir),
#                     save_float=True,
#                     save_8bit=False,
#                 )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#             hillshade_arr = self.get_hillshade(
#                 dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res, no_data=no_data
#             ).astype("float32")
#             if save_float:
#                 if (
#                     os.path.isfile(hillshade_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=hillshade_path,
#                         out_raster_arr=hillshade_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(hillshade_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     hillshade_8_bit_arr = self.float_to_8bit(
#                         float_arr=hillshade_arr,
#                         visualization=RVTVisualizationName.HILLSHADE,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=hillshade_8bit_path,
#                         out_raster_arr=hillshade_8_bit_arr,
#                         e_type=1,
#                     )
#             if save_shadow:
#                 shadow_arr = self.get_shadow(dem_arr=dem_arr, resolution=x_res)
#                 if (
#                     os.path.isfile(shadow_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=shadow_path,
#                         out_raster_arr=shadow_arr,
#                         no_data=np.nan,
#                     )
#             return True
#
#     def get_multi_hillshade(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution_x: float,
#         resolution_y: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         multi_hillshade_arr = rvt.vis.multi_hillshade(
#             dem=dem_arr,
#             resolution_x=resolution_x,
#             resolution_y=resolution_y,
#             number_of_directions=self.mhs_nr_dir,
#             sun_elevation=self.mhs_sun_el,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return multi_hillshade_arr
#
#     def save_multi_hillshade(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """
#         Calculates and saves Multidirectional hillshade from dem (dem_path) with default parameters.
#         If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255).
#         """
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.mhs_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.mhs_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_multi_hillshade: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_multi_hillshade: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             multi_hillshade_path = self.get_multi_hillshade_path(dem_path)
#             multi_hillshade_8bit_path = self.get_multi_hillshade_path(
#                 dem_path, is_8bit=True
#             )
#         else:
#             multi_hillshade_path = custom_dir / self.get_multi_hillshade_file_name(
#                 dem_path
#             )
#             multi_hillshade_8bit_path = custom_dir / self.get_multi_hillshade_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(multi_hillshade_8bit_path)
#                 and os.path.isfile(multi_hillshade_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(multi_hillshade_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(multi_hillshade_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#             if save_float:
#                 if (
#                     os.path.isfile(multi_hillshade_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     multi_hillshade_arr = self.get_multi_hillshade(
#                         dem_arr=dem_arr,
#                         resolution_x=x_res,
#                         resolution_y=y_res,
#                         no_data=no_data,
#                     ).astype("float32")
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=multi_hillshade_path,
#                         out_raster_arr=multi_hillshade_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(multi_hillshade_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     multi_hillshade_8bit_arr = self.float_to_8bit(
#                         float_arr=dem_arr,
#                         visualization=RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE,
#                         x_res=x_res,
#                         y_res=y_res,
#                         no_data=no_data,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=multi_hillshade_8bit_path,
#                         out_raster_arr=multi_hillshade_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_slrm(
#         self, dem_arr: npt.NDArray[Any], no_data: Optional[float] = None
#     ) -> npt.NDArray[Any]:
#         slrm_arr = rvt.vis.simple_local_relief_model(
#             dem=dem_arr,
#             radius_cell=self.slrm_rad_cell,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return slrm_arr
#
#     def save_slrm(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """Calculates and saves Simple local relief model from dem (dem_path) with default parameters.
#         If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255)."""
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.slrm_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.slrm_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_slrm: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_slrm: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             slrm_path = self.get_slrm_path(dem_path)
#             slrm_8bit_path = self.get_slrm_path(dem_path, is_8bit=True)
#         else:
#             slrm_path = custom_dir / self.get_slrm_file_name(dem_path)
#             slrm_8bit_path = custom_dir / self.get_slrm_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(slrm_8bit_path)
#                 and os.path.isfile(slrm_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(slrm_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(slrm_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.SIMPLE_LOCAL_RELIEF_MODEL,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             slrm_arr = self.get_slrm(dem_arr=dem_arr, no_data=no_data).astype("float32")
#             if save_float:
#                 if (
#                     os.path.isfile(slrm_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=slrm_path,
#                         out_raster_arr=slrm_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(slrm_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     slrm_8bit_arr = self.float_to_8bit(
#                         float_arr=slrm_arr,
#                         visualization=RVTVisualizationName.SIMPLE_LOCAL_RELIEF_MODEL,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=slrm_8bit_path,
#                         out_raster_arr=slrm_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_sky_view_factor(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution: float,
#         compute_svf: bool = True,
#         compute_asvf: bool = False,
#         compute_opns: bool = False,
#         no_data: Optional[float] = None,
#     ) -> Dict[str, npt.NDArray[Any]]:
#         dict_svf_asvf_opns = rvt.vis.horizon_visualizations(
#             dem=dem_arr,
#             resolution=resolution,
#             compute_svf=compute_svf,
#             compute_opns=compute_opns,
#             compute_asvf=compute_asvf,
#             number_of_directions=self.svf_n_dir,
#             maximum_search_radius=self.svf_r_max,
#             svf_noise=self.svf_noise,
#             direction_of_anisotropy=self.asvf_dir,
#             anisotropy_level=self.asvf_level,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return dict_svf_asvf_opns
#
#     def save_sky_view_factor(
#         self,
#         dem_path: Path,
#         save_svf: bool = True,
#         save_asvf: bool = False,
#         save_opns: bool = False,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """Calculates and saves Sky-view factor(save_svf=True), Anisotropic Sky-view factor(save_asvf=True) and
#         Positive Openness(save_opns=True) from dem (dem_path) with default parameters.
#         If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255)."""
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.svf_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.svf_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_sky_view_factor: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_sky_view_factor: dem_path doesn't exist!"
#             )
#
#         svf_path = Path()
#         asvf_path = Path()
#         opns_path = Path()
#         svf_8bit_path = Path()
#         asvf_8bit_path = Path()
#         opns_8bit_path = Path()
#         if custom_dir is None:
#             if save_svf:
#                 svf_path = self.get_svf_path(dem_path)
#                 svf_8bit_path = self.get_svf_path(dem_path, is_8bit=True)
#             if save_asvf:
#                 asvf_path = self.get_asvf_path(dem_path)
#                 asvf_8bit_path = self.get_asvf_path(dem_path, is_8bit=True)
#             if save_opns:
#                 opns_path = self.get_opns_path(dem_path)
#                 opns_8bit_path = self.get_opns_path(dem_path, is_8bit=True)
#         else:
#             if save_svf:
#                 svf_path = custom_dir / self.get_svf_file_name(dem_path)
#                 svf_8bit_path = custom_dir / self.get_svf_file_name(
#                     dem_path, is_8bit=True
#                 )
#             if save_asvf:
#                 asvf_path = custom_dir / self.get_asvf_file_name(dem_path)
#                 asvf_8bit_path = custom_dir / self.get_asvf_file_name(
#                     dem_path, is_8bit=True
#                 )
#             if save_opns:
#                 opns_path = custom_dir / self.get_opns_file_name(dem_path)
#                 opns_8bit_path = custom_dir / self.get_opns_file_name(
#                     dem_path, is_8bit=True
#                 )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(svf_path)
#                 and os.path.isfile(asvf_path)
#                 and os.path.isfile(opns_path)
#                 and os.path.isfile(svf_8bit_path)
#                 and os.path.isfile(asvf_8bit_path)
#                 and os.path.isfile(opns_8bit_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if (
#                 os.path.isfile(svf_path)
#                 and os.path.isfile(asvf_path)
#                 and os.path.isfile(opns_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif not save_float and save_8bit:
#             if (
#                 os.path.isfile(svf_8bit_path)
#                 and os.path.isfile(asvf_8bit_path)
#                 and os.path.isfile(opns_8bit_path)
#                 and not self.overwrite
#             ):
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             if save_svf:
#                 rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                     RVT_VISUALIZATION_NAME=RVTVisualizationName.SKY_VIEW_FACTOR,
#                     rvt_default=self,
#                     dem_path=Path(dem_path),
#                     output_dir_path=Path(custom_dir),
#                     save_float=save_float,
#                     save_8bit=save_8bit,
#                 )
#             if save_asvf:
#                 rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                     RVT_VISUALIZATION_NAME=RVTVisualizationName.ANISOTROPIC_SKY_VIEW_FACTOR,
#                     rvt_default=self,
#                     dem_path=Path(dem_path),
#                     output_dir_path=Path(custom_dir),
#                     save_float=save_float,
#                     save_8bit=save_8bit,
#                 )
#             if save_opns:
#                 rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                     RVT_VISUALIZATION_NAME=RVTVisualizationName.POSITIVE_OPENNESS,
#                     rvt_default=self,
#                     dem_path=Path(dem_path),
#                     output_dir_path=Path(custom_dir),
#                     save_float=save_float,
#                     save_8bit=save_8bit,
#                 )
#             return True
#         else:
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#             dict_svf_asvf_opns = self.get_sky_view_factor(
#                 dem_arr=dem_arr,
#                 resolution=x_res,
#                 compute_svf=save_svf,
#                 compute_asvf=save_asvf,
#                 compute_opns=save_opns,
#                 no_data=no_data,
#             )
#             if save_float:
#                 if save_svf:
#                     if (
#                         os.path.isfile(svf_path) and not self.overwrite
#                     ):  # file exists and overwrite=0
#                         pass
#                     else:  # svf_path, file doesn't exists or exists and overwrite=1
#                         save_raster(
#                             src_raster_path=dem_path,
#                             out_raster_path=svf_path,
#                             out_raster_arr=dict_svf_asvf_opns["svf"].astype("float32"),
#                             no_data=np.nan,
#                         )
#                 if save_asvf:
#                     if (
#                         os.path.isfile(asvf_path) and not self.overwrite
#                     ):  # file exists and overwrite=0
#                         pass
#                     else:  # asvf_path, file doesn't exists or exists and overwrite=1
#                         save_raster(
#                             src_raster_path=dem_path,
#                             out_raster_path=asvf_path,
#                             out_raster_arr=dict_svf_asvf_opns["asvf"].astype("float32"),
#                             no_data=np.nan,
#                         )
#                 if save_opns:
#                     if (
#                         os.path.isfile(opns_path) and not self.overwrite
#                     ):  # file exists and overwrite=0
#                         pass
#                     else:  # opns_path, file doesn't exists or exists and overwrite=1
#                         save_raster(
#                             src_raster_path=dem_path,
#                             out_raster_path=opns_path,
#                             out_raster_arr=dict_svf_asvf_opns["opns"].astype("float32"),
#                             no_data=np.nan,
#                         )
#             if save_8bit:
#                 if save_svf:
#                     if (
#                         os.path.isfile(svf_8bit_path) and not self.overwrite
#                     ):  # file exists and overwrite=0
#                         pass
#                     else:  # svf_8bit_path, file doesn't exists or exists and overwrite=1
#                         svf_8bit_arr = self.float_to_8bit(
#                             float_arr=dict_svf_asvf_opns["svf"],
#                             visualization=RVTVisualizationName.SKY_VIEW_FACTOR,
#                         )
#                         save_raster(
#                             src_raster_path=dem_path,
#                             out_raster_path=svf_8bit_path,
#                             out_raster_arr=svf_8bit_arr,
#                             e_type=1,
#                         )
#                 if save_asvf:
#                     if (
#                         os.path.isfile(asvf_8bit_path) and not self.overwrite
#                     ):  # file exists and overwrite=0
#                         pass
#                     else:  # asvf_8bit_path, file doesn't exists or exists and overwrite=1
#                         asvf_8bit_arr = self.float_to_8bit(
#                             float_arr=dict_svf_asvf_opns["asvf"],
#                             visualization=RVTVisualizationName.ANISOTROPIC_SKY_VIEW_FACTOR,
#                         )
#                         save_raster(
#                             src_raster_path=dem_path,
#                             out_raster_path=asvf_8bit_path,
#                             out_raster_arr=asvf_8bit_arr,
#                             e_type=1,
#                         )
#                 if save_opns:
#                     if (
#                         os.path.isfile(opns_8bit_path) and not self.overwrite
#                     ):  # file exists and overwrite=0
#                         pass
#                     else:  # opns_8bit_path, file doesn't exists or exists and overwrite=1
#                         opns_8bit_arr = self.float_to_8bit(
#                             float_arr=dict_svf_asvf_opns["opns"],
#                             visualization=RVTVisualizationName.POSITIVE_OPENNESS,
#                         )
#                         save_raster(
#                             src_raster_path=dem_path,
#                             out_raster_path=opns_8bit_path,
#                             out_raster_arr=opns_8bit_arr,
#                             e_type=1,
#                         )
#             return True
#
#     def get_neg_opns(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         dem_arr = -1 * dem_arr
#         dict_neg_opns = rvt.vis.horizon_visualizations(
#             dem=dem_arr,
#             resolution=resolution,
#             number_of_directions=self.svf_n_dir,
#             maximum_search_radius=self.svf_r_max,
#             svf_noise=self.svf_noise,
#             compute_svf=False,
#             compute_asvf=False,
#             compute_opns=True,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         neg_opns_arr = dict_neg_opns["opns"]
#         return neg_opns_arr
#
#     def save_neg_opns(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """
#         Calculates and saves Negative Openness from dem (dem_path) with default parameters. If custom_dir is None
#         it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255).
#         """
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.neg_opns_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.neg_opns_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_neg_opns: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_neg_opns: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             neg_opns_path = self.get_neg_opns_path(dem_path)
#             neg_opns_8bit_path = self.get_neg_opns_path(dem_path, is_8bit=True)
#         else:
#             neg_opns_path = custom_dir / self.get_neg_opns_file_name(dem_path)
#             neg_opns_8bit_path = custom_dir / self.get_neg_opns_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(neg_opns_8bit_path)
#                 and os.path.isfile(neg_opns_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(neg_opns_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(neg_opns_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.NEGATIVE_OPENNESS,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#
#             neg_opns_arr = self.get_neg_opns(
#                 dem_arr=dem_arr, resolution=x_res, no_data=no_data
#             ).astype("float32")
#             if save_float:
#                 if (
#                     os.path.isfile(neg_opns_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=neg_opns_path,
#                         out_raster_arr=neg_opns_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(neg_opns_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     neg_opns_8bit_arr = self.float_to_8bit(
#                         float_arr=neg_opns_arr,
#                         visualization=RVTVisualizationName.NEGATIVE_OPENNESS,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=neg_opns_8bit_path,
#                         out_raster_arr=neg_opns_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_sky_illumination(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         sky_illumination_arr = rvt.vis.sky_illumination(
#             dem=dem_arr,
#             resolution=resolution,
#             sky_model=self.sim_sky_mod,
#             compute_shadow=bool(self.sim_compute_shadow),
#             max_fine_radius=self.sim_shadow_dist,
#             num_directions=self.sim_nr_dir,
#             shadow_az=self.sim_shadow_az,
#             shadow_el=self.sim_shadow_el,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return sky_illumination_arr  # type: ignore
#
#     def save_sky_illumination(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """Calculates and saves Sky illumination from dem (dem_path) with default parameters. If custom_dir is None
#         it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255)."""
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.sim_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.sim_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_sky_illumination: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_sky_illumination: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             sky_illumination_path = self.get_sky_illumination_path(dem_path)
#             sky_illumination_8bit_path = self.get_sky_illumination_path(
#                 dem_path, is_8bit=True
#             )
#         else:
#             sky_illumination_path = custom_dir / self.get_sky_illumination_file_name(
#                 dem_path
#             )
#             sky_illumination_8bit_path = (
#                 custom_dir / self.get_sky_illumination_file_name(dem_path, is_8bit=True)
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(sky_illumination_8bit_path)
#                 and os.path.isfile(sky_illumination_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(sky_illumination_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(sky_illumination_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.SKY_ILLUMINATION,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#
#             sky_illumination_arr = self.get_sky_illumination(
#                 dem_arr=dem_arr, resolution=x_res, no_data=no_data
#             ).astype("float32")
#             if save_float:
#                 if (
#                     os.path.isfile(sky_illumination_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=sky_illumination_path,
#                         out_raster_arr=sky_illumination_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(sky_illumination_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     sky_illumination_8bit_arr = self.float_to_8bit(
#                         float_arr=sky_illumination_arr,
#                         visualization=RVTVisualizationName.SKY_ILLUMINATION,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=sky_illumination_8bit_path,
#                         out_raster_arr=sky_illumination_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_local_dominance(
#         self, dem_arr: npt.NDArray[Any], no_data: Optional[float] = None
#     ) -> npt.NDArray[Any]:
#         local_dominance_arr = rvt.vis.local_dominance(
#             dem=dem_arr,
#             minimum_radius=self.ld_min_rad,
#             maximum_radius=self.ld_max_rad,
#             radial_distance_step=self.ld_rad_inc,
#             angular_step=self.ld_anglr_res,
#             observer_height=self.ld_observer_h,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return local_dominance_arr
#
#     def save_local_dominance(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """
#         Calculates and saves Local dominance from dem (dem_path) with default parameters. If custom_dir is None
#         it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255).
#         """
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.ld_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.ld_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_local_dominance: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_local_dominance: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             local_dominance_path = self.get_local_dominance_path(dem_path)
#             local_dominance_8bit_path = self.get_local_dominance_path(
#                 dem_path, is_8bit=True
#             )
#         else:
#             local_dominance_path = custom_dir / self.get_local_dominance_file_name(
#                 dem_path
#             )
#             local_dominance_8bit_path = custom_dir / self.get_local_dominance_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(local_dominance_8bit_path)
#                 and os.path.isfile(local_dominance_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(local_dominance_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(local_dominance_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.LOCAL_DOMINANCE,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             local_dominance_arr = self.get_local_dominance(
#                 dem_arr=dem_arr, no_data=no_data
#             ).astype("float32")
#             if save_float:
#                 if (
#                     os.path.isfile(local_dominance_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=local_dominance_path,
#                         out_raster_arr=local_dominance_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(local_dominance_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     local_dominance_8bit_arr = self.float_to_8bit(
#                         float_arr=local_dominance_arr,
#                         visualization=RVTVisualizationName.LOCAL_DOMINANCE,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=local_dominance_8bit_path,
#                         out_raster_arr=local_dominance_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_msrm(
#         self,
#         dem_arr: npt.NDArray[Any],
#         resolution: float,
#         no_data: Optional[float] = None,
#     ) -> npt.NDArray[Any]:
#         msrm_arr = rvt.vis.multi_scale_relief_model(
#             dem=dem_arr,
#             resolution=resolution,
#             minimum_feature_size=self.msrm_feature_min,
#             maximum_feature_size=self.msrm_feature_max,
#             scaling_factor=self.msrm_scaling_factor,
#             vertical_exaggeration_factor=self.vertical_exaggeration_factor,
#             no_data=no_data,
#         )
#         return msrm_arr
#
#     def save_msrm(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """
#         Calculates and saves Multi-scale relief model from dem (dem_path) with default parameters.
#         If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
#         if save_8bit is True method creates GTiff with bytescaled values (0-255).
#         """
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.msrm_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.msrm_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_msrm: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_msrm: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             msrm_path = self.get_msrm_path(dem_path)
#             msrm_8bit_path = self.get_msrm_path(dem_path, is_8bit=True)
#         else:
#             msrm_path = custom_dir / self.get_msrm_file_name(dem_path)
#             msrm_8bit_path = custom_dir / self.get_msrm_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(msrm_8bit_path)
#                 and os.path.isfile(msrm_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(msrm_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(msrm_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit:  # tile by tile
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.MULTI_SCALE_RELIEF_MODEL,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#             x_res = dict_arr_res["resolution"][0]
#             y_res = dict_arr_res["resolution"][1]
#
#             msrm_arr = self.get_msrm(
#                 dem_arr=dem_arr, resolution=x_res, no_data=no_data
#             ).astype("float32")
#             if save_float:
#                 if (
#                     os.path.isfile(msrm_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=msrm_path,
#                         out_raster_arr=msrm_arr,
#                         no_data=np.nan,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(msrm_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     msrm_8bit_arr = self.float_to_8bit(
#                         float_arr=msrm_arr,
#                         visualization=RVTVisualizationName.MULTI_SCALE_RELIEF_MODEL,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=msrm_8bit_path,
#                         out_raster_arr=msrm_8bit_arr,
#                         e_type=1,
#                     )
#             return True
#
#     def get_mstp(
#         self, dem_arr: npt.NDArray[Any], no_data: Optional[float] = None
#     ) -> npt.NDArray[Any]:
#         mstp_arr = rvt.vis.multi_scale_topographic_position(
#             dem=dem_arr,
#             local_scale=self.mstp_local_scale,
#             meso_scale=self.mstp_meso_scale,
#             broad_scale=self.mstp_broad_scale,
#             lightness=self.mstp_lightness,
#             no_data=no_data,
#         )
#         return mstp_arr
#
#     def save_mstp(
#         self,
#         dem_path: Path,
#         custom_dir: Optional[Path] = None,
#         save_float: Optional[bool] = None,
#         save_8bit: Optional[bool] = None,
#     ) -> bool:
#         """Calculates and saves Multi-scale topographic position from dem (dem_path) with default parameters.
#         If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
#         overwrite file (overwrite=1) or not (overwrite=0)."""
#
#         # if save_float is None it takes boolean from default (self)
#         if save_float is None:
#             save_float = bool(self.mstp_save_float)
#         # if save_8bit is None it takes boolean from default (self)
#         if save_8bit is None:
#             save_8bit = bool(self.mstp_save_8bit)
#
#         if not save_float and not save_8bit:
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_mstp: Both save_float and save_8bit are False,"
#                 " at least one of them has to be True!"
#             )
#
#         if not os.path.isfile(dem_path):
#             raise Exception(
#                 "rvt.default.RVTVisualizationFactory.save_mstp: dem_path doesn't exist!"
#             )
#
#         if custom_dir is None:
#             mstp_path = self.get_mstp_path(dem_path)
#             mstp_8bit_path = self.get_mstp_path(dem_path, is_8bit=True)
#         else:
#             mstp_path = custom_dir / self.get_mstp_file_name(dem_path)
#             mstp_8bit_path = custom_dir / self.get_mstp_file_name(
#                 dem_path, is_8bit=True
#             )
#
#         # if file already exists and overwrite=0
#         if save_float and save_8bit:
#             if (
#                 os.path.isfile(mstp_8bit_path)
#                 and os.path.isfile(mstp_path)
#                 and not self.overwrite
#             ):
#                 return False
#         elif save_float and not save_8bit:
#             if os.path.isfile(mstp_path) and not self.overwrite:
#                 return False
#         elif not save_float and not save_8bit:
#             if os.path.isfile(mstp_8bit_path) and not self.overwrite:
#                 return False
#
#         dem_size = get_raster_size(raster_path=dem_path)
#         if (
#             dem_size[0] * dem_size[1] > self.compute_multiple_tile_limit
#         ):  # tile by tile calculation
#             if custom_dir is None:
#                 custom_dir = Path(dem_path).parent
#             rvt.tile.save_RVT_VISUALIZATION_NAME_tile_by_tile(
#                 RVT_VISUALIZATION_NAME=RVTVisualizationName.MULTI_SCALE_TOPOGRAPHIC_POSITION,
#                 rvt_default=self,
#                 dem_path=Path(dem_path),
#                 output_dir_path=Path(custom_dir),
#                 save_float=save_float,
#                 save_8bit=save_8bit,
#             )
#             return True
#         else:  # singleprocess
#             dict_arr_res = get_raster_arr(raster_path=dem_path)
#             dem_arr = dict_arr_res["array"]
#             no_data = dict_arr_res["no_data"]
#
#             mstp_arr = self.get_mstp(dem_arr=dem_arr, no_data=no_data)
#
#             if save_float:
#                 if (
#                     os.path.isfile(mstp_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=mstp_path,
#                         out_raster_arr=mstp_arr,
#                         no_data=np.nan,
#                         e_type=6,
#                     )
#             if save_8bit:
#                 if (
#                     os.path.isfile(mstp_8bit_path) and not self.overwrite
#                 ):  # file exists and overwrite=0
#                     pass
#                 else:
#                     mstp_8bit_arr = self.float_to_8bit(
#                         float_arr=mstp_arr,
#                         visualization=RVTVisualizationName.MULTI_SCALE_TOPOGRAPHIC_POSITION,
#                     )
#                     save_raster(
#                         src_raster_path=dem_path,
#                         out_raster_path=mstp_8bit_path,
#                         out_raster_arr=mstp_8bit_arr,
#                         no_data=np.nan,
#                         e_type=1,
#                     )
#
#             return True
#
#     def save_visualizations(
#         self, dem_path: Path, custom_dir: Optional[Path] = None
#     ) -> None:
#         """Save all visualizations where self.'visualization'_compute = True also saves float where self.'visualization'
#         _save_float = True and 8bit where self.'visualization'_save_8bit = True. In the end method creates log file."""
#         start_time = time.time()
#         if self.slp_compute:
#             self.save_slope(dem_path, custom_dir=custom_dir)
#         if self.hs_compute:
#             self.save_hillshade(dem_path, custom_dir=custom_dir)
#         if self.mhs_compute:
#             self.save_multi_hillshade(dem_path, custom_dir=custom_dir)
#         if self.slrm_compute:
#             self.save_slrm(dem_path, custom_dir=custom_dir)
#         if self.svf_compute or self.asvf_compute or self.pos_opns_compute:
#             self.save_sky_view_factor(
#                 dem_path,
#                 save_svf=bool(self.svf_compute),
#                 save_asvf=bool(self.asvf_compute),
#                 save_opns=bool(self.pos_opns_compute),
#                 custom_dir=custom_dir,
#             )
#         if self.neg_opns_compute:
#             self.save_neg_opns(dem_path, custom_dir=custom_dir)
#         if self.sim_compute:
#             self.save_sky_illumination(dem_path, custom_dir=custom_dir)
#         if self.ld_compute:
#             self.save_local_dominance(dem_path, custom_dir=custom_dir)
#         if self.msrm_compute:
#             self.save_msrm(dem_path, custom_dir=custom_dir)
#         if self.mstp_compute:
#             self.save_mstp(dem_path, custom_dir=custom_dir)
#         end_time = time.time()
#         compute_time = end_time - start_time
#         self.create_log_file(
#             dem_path=dem_path, custom_dir=custom_dir, compute_time=compute_time
#         )
#
#     def calculate_visualization(
#         self,
#         visualization: RVTVisualizationName,
#         dem: npt.NDArray[Any],
#         resolution_x: float,
#         resolution_y: float,
#         no_data: Optional[float] = None,
#         save_float: bool = True,
#         save_8bit: bool = False,
#     ) -> Tuple[
#         Optional[npt.NDArray], Optional[npt.NDArray]
#     ]:  # tuple[vis_float_arr, vis_8bit_arr]
#         vis_arr = None
#         vis_float_arr = None
#         vis_8bit_arr = None
#         if visualization == RVTVisualizationName.SLOPE:
#             vis_arr = self.get_slope(
#                 dem_arr=dem,
#                 resolution_x=resolution_x,
#                 resolution_y=resolution_y,
#                 no_data=no_data,
#             )
#         elif visualization == RVTVisualizationName.SHADOW:
#             vis_arr = self.get_shadow(
#                 dem_arr=dem, resolution=resolution_x, no_data=no_data
#             )
#         elif visualization == RVTVisualizationName.HILLSHADE:
#             vis_arr = self.get_hillshade(
#                 dem_arr=dem,
#                 resolution_x=resolution_x,
#                 resolution_y=resolution_y,
#                 no_data=no_data,
#             )
#         elif visualization == RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE:
#             vis_arr = self.get_multi_hillshade(
#                 dem_arr=dem,
#                 resolution_x=resolution_x,
#                 resolution_y=resolution_y,
#                 no_data=no_data,
#             )
#         elif visualization == RVTVisualizationName.SIMPLE_LOCAL_RELIEF_MODEL:
#             vis_arr = self.get_slrm(dem_arr=dem, no_data=no_data)
#         elif visualization == RVTVisualizationName.SKY_VIEW_FACTOR:
#             vis_arr = self.get_sky_view_factor(
#                 dem_arr=dem,
#                 resolution=resolution_x,
#                 compute_svf=True,
#                 compute_asvf=False,
#                 compute_opns=False,
#                 no_data=no_data,
#             )["svf"]
#         elif visualization == RVTVisualizationName.ANISOTROPIC_SKY_VIEW_FACTOR:
#             vis_arr = self.get_sky_view_factor(
#                 dem_arr=dem,
#                 resolution=resolution_x,
#                 compute_svf=False,
#                 compute_asvf=True,
#                 compute_opns=False,
#                 no_data=no_data,
#             )["asvf"]
#         elif visualization == RVTVisualizationName.POSITIVE_OPENNESS:
#             vis_arr = self.get_sky_view_factor(
#                 dem_arr=dem,
#                 resolution=resolution_x,
#                 compute_svf=False,
#                 compute_asvf=False,
#                 compute_opns=True,
#                 no_data=no_data,
#             )["opns"]
#         elif visualization == RVTVisualizationName.NEGATIVE_OPENNESS:
#             vis_arr = self.get_neg_opns(
#                 dem_arr=dem, resolution=resolution_x, no_data=no_data
#             )
#         elif visualization == RVTVisualizationName.SKY_ILLUMINATION:
#             vis_arr = self.get_sky_illumination(
#                 dem_arr=dem, resolution=resolution_x, no_data=no_data
#             )
#         elif visualization == RVTVisualizationName.LOCAL_DOMINANCE:
#             vis_arr = self.get_local_dominance(dem_arr=dem, no_data=no_data)
#         elif visualization == RVTVisualizationName.MULTI_SCALE_RELIEF_MODEL:
#             vis_arr = self.get_msrm(
#                 dem_arr=dem, resolution=resolution_x, no_data=no_data
#             )
#         elif visualization == RVTVisualizationName.MULTI_SCALE_TOPOGRAPHIC_POSITION:
#             vis_arr = self.get_mstp(dem_arr=dem, no_data=no_data)
#         if save_float:
#             vis_float_arr = vis_arr
#         if save_8bit:
#             if visualization == RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE:
#                 vis_8bit_arr = self.float_to_8bit(
#                     float_arr=dem,
#                     visualization=visualization,
#                     x_res=resolution_x,
#                     y_res=resolution_y,
#                     no_data=no_data,
#                 )
#             else:
#                 assert vis_arr is not None
#                 vis_8bit_arr = self.float_to_8bit(
#                     float_arr=vis_arr, visualization=visualization
#                 )
#         return vis_float_arr, vis_8bit_arr
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
#
# def get_raster_size(raster_path: Path, band_number: int = 1) -> Tuple[float, float]:
#     # TODO: Refactor and replace with rasterio
#     """Opens raster path and returns selected band size.
#
#     Parameters
#     ----------
#     raster_path : str
#         Path to raster.
#     band_number : int
#         Selected band number.
#
#     Returns
#     -------
#     tuple(x_size, y_size)
#     """
#     data_set = gdal.Open(raster_path)  # Open dem raster
#     band = data_set.GetRasterBand(band_number)
#     x_size = band.XSize  # number of columns
#     y_size = band.YSize  # number of rows
#     del band
#     data_set = None  # close data_set
#     return x_size, y_size
#
#
# def save_raster(
#     src_raster_path: Path,
#     out_raster_path: Path,
#     out_raster_arr: npt.NDArray[Any],
#     no_data: Optional[float] = None,
#     e_type: int = 6,
# ) -> None:
#     # TODO: Refactor and replace with rasterio
#     """Saves raster array (out_rast_arr) to out_raster_path (GTiff), using src_rast_path information.
#
#     Parameters
#     ----------
#     src_raster_path : str
#         Path to source raster.
#     out_raster_path : str
#         Path to new file, where to save raster (GTiff).
#     out_raster_arr : np.array (2D - one band_number, 3D - multiple bands)
#         Array with raster data.
#     no_data : float
#         Value that represents no data pixels.
#     e_type : GDALDataType
#         https://gdal.org/api/raster_c_api.html#_CPPv412GDALDataType, (GDT_Float32 = 6, GDT_UInt8 = 1, ...)
#     """
#     src_data_set = gdal.Open(src_raster_path)
#     gtiff_driver = gdal.GetDriverByName("GTiff")
#     if len(out_raster_arr.shape) == 2:  # 2D array, one band_number
#         out_data_set = gtiff_driver.Create(
#             out_raster_path,
#             xsize=out_raster_arr.shape[1],
#             ysize=out_raster_arr.shape[0],
#             bands=1,
#             eType=e_type,  # eType: 6 = GDT_Float32
#             options=["COMPRESS=LZW"],
#         )
#         out_data_set.SetProjection(src_data_set.GetProjection())
#         out_data_set.SetGeoTransform(src_data_set.GetGeoTransform())
#         out_data_set.GetRasterBand(1).WriteArray(out_raster_arr)
#         if no_data is not None:
#             out_data_set.GetRasterBand(1).SetNoDataValue(no_data)
#
#     elif len(out_raster_arr.shape) == 3:  # 3D array, more bands
#         out_data_set = gtiff_driver.Create(
#             out_raster_path,
#             xsize=out_raster_arr.shape[2],
#             ysize=out_raster_arr.shape[1],
#             bands=out_raster_arr.shape[0],
#             eType=e_type,  # eType: 6 = GDT_Float32
#             options=["COMPRESS=LZW"],
#         )
#         out_data_set.SetProjection(src_data_set.GetProjection())
#         out_data_set.SetGeoTransform(src_data_set.GetGeoTransform())
#         for i_band in range(out_raster_arr.shape[0]):
#             out_data_set.GetRasterBand(i_band + 1).WriteArray(
#                 out_raster_arr[i_band, :, :]
#             )
#         if no_data is not None:
#             out_data_set.GetRasterBand(1).SetNoDataValue(no_data)
#     else:
#         raise Exception(
#             "rvt.default.save_raster: You have to input 2D or 3D numpy array!"
#         )
#     out_data_set.FlushCache()
#     src_data_set = None  # Close source data set
#     out_data_set = None  # Close output data set
