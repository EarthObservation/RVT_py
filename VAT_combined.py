"""
This python script is for calculating VAT Combined blender combination
(1st layer: VAT general (with opacity set in parameters, blend mode normal), 2nd layer: VAT flat)
on multiple DEMs with multiple processes.
Parameters are:
    input_dir_path (dir that contains GeoTIFFs of DEMs)
    output_dir_path (dir where to store calculated GeoTIFFS of VAT combined visualizations)
    general_opacity (opacity of VAT general)
    vat_combination_json_path (path to VAT combination JSON file)
    terrains_sett_json_path (path to terrains settings JSON file)
    nr_processes (number of parallel processes which are calculating and saving output VAT combined)
    save_float (if saves float)
    save_8bit (if saves 8bit)
    save_VAT_general (if saves VAT general)
    save_VAT_flat (if saves VAT flat)
"""
import rvt.visualizations
import rvt.blender
import rvt.default
import os
import multiprocessing as mp

input_dir_path = "input_dir"
output_dir_path = "output_dir"
general_opacity = 50
vat_combination_json_path = "settings/blender_VAT.json"
terrains_sett_json_path = "settings/default_terrains_settings.json"
nr_processes = 2
save_float = True
save_8bit = True
save_VAT_general = True
save_VAT_flat = True


def combined_VAT(
    input_dir_path,
    output_dir_path,
    general_opacity,
    vat_combination_json_path=None,
    terrains_sett_json_path=None,
    nr_processes=7,
    save_float=True,
    save_8bit=False,
    save_VAT_general=False,
    save_VAT_flat=False,
):
    if not save_float and not save_8bit:
        raise Exception("save_float and save_8bit are both False!")

    if (
        vat_combination_json_path is None
    ):  # Če ni podan path do vat_comb se smatra da je v settings
        vat_combination_json_path = os.path.abspath(
            os.path.join("settings", "blender_VAT.json")
        )

    if (
        terrains_sett_json_path is None
    ):  # Če ni podan path do terrain sett se smatra da je v settings
        terrains_sett_json_path = os.path.abspath(
            os.path.join("settings", "default_terrains_settings.json")
        )

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
    vat_combination_1.read_from_json_file(vat_combination_json_path)
    vat_combination_2.read_from_json_file(vat_combination_json_path)

    # create terrains settings and read it from JSON
    terrains_settings = rvt.blend.TerrainsSettings()
    terrains_settings.read_parameters_from_file(terrains_sett_json_path)

    # select single terrain settings for general and flat
    terrain_1 = terrains_settings.select_terrain_settings_by_name(
        "general"
    )  # VAT general
    terrain_2 = terrains_settings.select_terrain_settings_by_name("flat")  # VAT flat

    # apply terrain settings on default and combination
    terrain_1.apply_terrain(
        default=default_1, combination=vat_combination_1
    )  # VAT general
    terrain_2.apply_terrain(
        default=default_2, combination=vat_combination_2
    )  # VAT flat

    # prepare for multiprocessing
    dem_list = os.listdir(input_dir_path)
    input_process_list = []
    for input_dem_name in dem_list:
        if ".tif" not in input_dem_name:  # preskoči če se file ne konča .tif
            continue
        input_dem_path = os.path.join(input_dir_path, input_dem_name)
        out_name = "{}_Archaeological_(VAT_combined)_opac{}.tif".format(
            input_dem_name.rstrip(".tif"), general_opacity
        )
        out_comb_vat_path = os.path.abspath(os.path.join(output_dir_path, out_name))
        out_comb_vat_8bit_path = out_comb_vat_path.rstrip(".tif") + "_8bit.tif"
        if (
            save_8bit
            and os.path.isfile(out_comb_vat_8bit_path)
            and save_float
            and os.path.isfile(out_comb_vat_path)
        ):
            print("{} already exists!".format(out_comb_vat_path))
            print("{} already exists!".format(out_comb_vat_8bit_path))
            continue
        elif (
            save_float and os.path.isfile(out_comb_vat_path) and not save_8bit
        ):  # output VAT comb already exists
            print("{} already exists!".format(out_comb_vat_path))
            continue
        elif save_8bit and os.path.isfile(out_comb_vat_8bit_path) and not save_float:
            print("{} already exists!".format(out_comb_vat_8bit_path))
            continue
        general_combination = vat_combination_1
        flat_combination = vat_combination_2
        general_default = default_1
        flat_default = default_2
        out_comb_vat_general_name = "{}_Archaeological_(VAT_general).tif".format(
            input_dem_name.rstrip(".tif")
        )
        out_comb_vat_general_path = os.path.abspath(
            os.path.join(output_dir_path, out_comb_vat_general_name)
        )
        out_comb_vat_flat_name = "{}_Archaeological_(VAT_flat).tif".format(
            input_dem_name.rstrip(".tif")
        )
        out_comb_vat_flat_path = os.path.abspath(
            os.path.join(output_dir_path, out_comb_vat_flat_name)
        )
        input_process_list.append(
            (
                general_combination,
                flat_combination,
                general_default,
                flat_default,
                input_dem_path,
                out_comb_vat_path,
                general_opacity,
                save_float,
                save_8bit,
                save_VAT_general,
                out_comb_vat_general_path,
                save_VAT_flat,
                out_comb_vat_flat_path,
            )
        )

    # multiprocessing
    with mp.Pool(nr_processes) as p:
        realist = [
            p.apply_async(compute_save_VAT_combined, r) for r in input_process_list
        ]
        for result in realist:
            print(result.get())


# function which is multiprocessing
def compute_save_VAT_combined(
    general_combination,
    flat_combination,
    general_default,
    flat_default,
    input_dem_path,
    out_comb_vat_path,
    general_transparency,
    save_float,
    save_8bit,
    save_VAT_general,
    out_comb_vat_general_path,
    save_VAT_flat,
    out_comb_vat_flat_path,
):
    dict_arr_res_nd = rvt.default.get_raster_arr(raster_path=input_dem_path)

    # create and blend VAT general
    general_combination.add_dem_arr(
        dem_arr=dict_arr_res_nd["array"],
        dem_resolution=dict_arr_res_nd["resolution"][0],
    )
    if save_VAT_general:
        general_combination.add_dem_path(dem_path=input_dem_path)
        vat_arr_1 = general_combination.render_all_images(
            save_render_path=out_comb_vat_general_path,
            save_float=save_float,
            save_8bit=save_8bit,
            default=general_default,
            no_data=dict_arr_res_nd["no_data"],
        )
    else:
        vat_arr_1 = general_combination.render_all_images(
            default=general_default, no_data=dict_arr_res_nd["no_data"]
        )

    # create and blend VAT flat
    flat_combination.add_dem_arr(
        dem_arr=dict_arr_res_nd["array"],
        dem_resolution=dict_arr_res_nd["resolution"][0],
    )
    if save_VAT_flat:
        flat_combination.add_dem_path(dem_path=input_dem_path)
        vat_arr_2 = flat_combination.render_all_images(
            save_render_path=out_comb_vat_flat_path,
            save_float=save_float,
            save_8bit=save_8bit,
            default=flat_default,
            no_data=dict_arr_res_nd["no_data"],
        )
    else:
        vat_arr_2 = flat_combination.render_all_images(
            default=flat_default, no_data=dict_arr_res_nd["no_data"]
        )

    # create combination which blends VAT flat and VAT general together
    combination = rvt.blend.BlenderCombination()
    combination.create_layer(
        vis_method="VAT general",
        image=vat_arr_1,
        normalization="Value",
        minimum=0,
        maximum=1,
        blend_mode="Normal",
        opacity=general_transparency,
    )
    combination.create_layer(
        vis_method="VAT flat",
        image=vat_arr_2,
        normalization="Value",
        minimum=0,
        maximum=1,
        blend_mode="Normal",
        opacity=100,
    )

    combination.add_dem_path(dem_path=input_dem_path)

    # VAT combined
    combination.render_all_images(
        save_render_path=out_comb_vat_path,
        save_visualizations=False,
        save_float=save_float,
        save_8bit=save_8bit,
        no_data=dict_arr_res_nd["no_data"],
    )
    out_comb_vat_8bit_path = out_comb_vat_path.rstrip("tif") + "_8bit.tif"
    if save_float and save_8bit:
        return "{} and {} successfully calculated and saved!".format(
            out_comb_vat_path, out_comb_vat_8bit_path
        )
    elif save_float:
        return "{} successfully calculated and saved!".format(out_comb_vat_path)
    elif save_8bit:
        return "{} successfully calculated and saved!".format(out_comb_vat_8bit_path)


# Program start
if __name__ == "__main__":
    combined_VAT(
        input_dir_path=input_dir_path,
        output_dir_path=output_dir_path,
        general_opacity=general_opacity,
        vat_combination_json_path=vat_combination_json_path,
        terrains_sett_json_path=terrains_sett_json_path,
        nr_processes=nr_processes,
        save_float=save_float,
        save_8bit=save_8bit,
        save_VAT_general=save_VAT_general,
        save_VAT_flat=save_VAT_flat,
    )
