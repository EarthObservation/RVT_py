"""
GUI for tiled processing with RVT_py
Created on 6 May 2024
@author: Nejc Čož, ZRC SAZU, Novi trg 2, 1000 Ljubljana, Slovenia
"""

from pathlib import Path
from tkinter import Tk, filedialog

import geopandas as gpd
import ipywidgets as widgets
import traitlets
from IPython.display import display

import grid_tools as gt
from tiled_multiprocess import tiled_blending


class SelectFilesButton(widgets.Button):
    def __init__(self):
        super(SelectFilesButton, self).__init__()

        # Add traits
        self.add_traits(files=traitlets.traitlets.List())

        # Button options
        self.description = "Select file"
        self.style.button_color = None

        # Set on click behavior.
        self.on_click(self.select_files)

    @staticmethod
    def select_files(b):
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)

        selected_files = filedialog.askopenfilename(
            title="Select input files",
            filetypes=[("GeoTIF", "*.tif;*.tiff"), ("VRT", "*.vrt")],
            multiple=True
        )

        n_files = len(selected_files)
        # IF THERE WAS AT LEAST ONE FILE SELECTED
        if n_files > 0 and any(element != "" for element in selected_files):
            # Change button style
            b.description = f"{n_files} File selected" if n_files == 1 else f"{n_files} Files selected"
            b.style.button_color = "lightgreen"

            b.files = selected_files


# Button - opens dialog window (select file/s)
b_file_select = SelectFilesButton()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
dropdown_options = widgets.SelectMultiple(
    options=[
        ("VAT combined - 8 bit", 'vat_combined_8bit'),
        # "VAT_combined_3B",
        ('SLRM', 'SLRM'),
        ('e2MSTP', 'e2MSTP'),
        # ('e3MSTP', 'e3MSTP'),
        ('e4MSTP', 'e4MSTP')
    ],
    value=[],
    rows=10,
    disabled=False
)

# RUN BUTTON ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
button_run_adaf = widgets.Button(
    description="Run RVT",
    # layout={'width': 'auto', 'border': '1px solid black'},
    tooltip='Description'
)


def on_button_clicked(b):
    output_widget.clear_output()
    button_run_adaf.disabled = True

    list_tifs = b_file_select.files

    blend_types = list(dropdown_options.value)

    # CONSTANT INPUT PARAMETERS
    tile_size = 1000  # in PIXELS!
    refgrid_existing = False  # Change to true, to avoid re-creating refgrid
    nr_processes = 4  # How many CPUs to use?

    with output_widget:
        print(list_tifs)
        print(blend_types)

        for in_file in list_tifs:
            in_file = Path(in_file)
            # (1b) Instead of folder, the tif/vrt is given
            input_vrt = in_file.as_posix()
            ds = in_file.parent

            print("Start --- " + ds.name)

            # (2) and (3) if refgrid doesn't exist, else read it from file
            if refgrid_existing:
                # I have added this option for Noise Mapping, because vrt and refgrid are filtered to contain
                # only the tiles with archeology.
                print("* refgrid exists; Reading from file...")
                refgrid_name = list(ds.glob("*_refgrid*.gpkg"))[0]
                tiles_extents = gpd.read_file(refgrid_name)
            else:
                # (2) To filter we need polygon covering valid data
                vdm_file = list(ds.glob("*_validDataMask*"))

                if vdm_file:
                    print("* validDataMask exists; Removing and creating new...")
                    Path(vdm_file[0]).unlink()
                else:
                    # If it doesn't exist, try creating from VRT
                    print("* validDataMask doesn't exists; Creating...")

                valid_data_outline = gt.poly_from_valid(input_vrt, save_gpkg=True)

                # (3) Create reference grid, filter it and save it to disk
                print("Create REF GRID")

                refgrid_name = input_vrt[:-4] + "_refgrid.gpkg"
                if Path(refgrid_name).exists():
                    Path(refgrid_name).unlink()

                tiles_extents = gt.bounding_grid(input_vrt, tile_size, tag=True)
                tiles_extents = gt.filter_by_outline(
                    tiles_extents, valid_data_outline.as_posix(),
                    save_gpkg=True, save_path=refgrid_name
                )

            # Run tiled blending (visualizations calculated on the go, stored in memory)
            tiles_list = tiles_extents[["minx", "miny", "maxx", "maxy"]].values.tolist()

            tiled_blending(
                blend_types=blend_types,
                input_vrt_path=in_file,
                tiles_list=tiles_list,
                nr_processes=nr_processes
            )

    button_run_adaf.disabled = False


button_run_adaf.on_click(on_button_clicked)

# DISPLAY OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
output_widget = widgets.Output()

# DISPLAY WIDGET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
display(
    widgets.VBox(
        [
            widgets.HTML(value=f"<b>RVT - TILED PROCESSING FOR LARGE FILES</b>"),
            b_file_select,
            widgets.Label("Select visualizations:"),
            dropdown_options,
            button_run_adaf,
            output_widget
        ]
    )
)
