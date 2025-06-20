"""
GUI for tiled processing with RVT_py
Created on 6 May 2024
@author: Nejc Čož, ZRC SAZU, Novi trg 2, 1000 Ljubljana, Slovenia
"""
from tkinter import Tk, filedialog

import ipywidgets as widgets
import traitlets
from IPython.display import display

from tiled_multiprocess import run_main


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
dropdown_options_blends = widgets.SelectMultiple(
    options=[
        ("VAT combined", 'vat_combined'),
        # "VAT_combined_3B",
        ('e2MSTP', 'e2MSTP'),
        # ('e3MSTP', 'e3MSTP'),
        ('e4MSTP', 'e4MSTP'),
        ('RRIM', 'rrim'),
        ('CRIM', 'crim')
    ],
    value=[],
    rows=10,
    disabled=False
)

dropdown_options_visualizations = widgets.SelectMultiple(
    options=[
        ('Slope', 'slp'),
        ('Hillshade', 'hs'),
        # ('multi_hillshade', 'multi_hillshade'),
        ('SLRM', 'slrm'),
        ('Sky view factor', 'svf'),
        ('Openness', 'opns'),
        ('Negative openness', 'neg_opns'),
        ('Local dominance', 'ld'),
        # ('sky_illumination', 'sky_illumination'),
        # ('shadow_horizon', 'shadow_horizon'),
        # ('MSRM', 'msrm'),
        ('MSTP', 'mstp')
    ],
    value=[],
    rows=10,
    disabled=False
)

# Convert to 8bit image
save_float_checkbox = widgets.Checkbox(
    value=False,
    description='Save float',
    disabled=False,
    indent=False
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

    vis_types = list(dropdown_options_visualizations.value)
    blend_types = list(dropdown_options_blends.value)

    save_float = save_float_checkbox.value

    with output_widget:
        print(list_tifs)
        print(vis_types)
        print(blend_types)
        # print(str(save_float_checkbox.value))

        run_main(list_tifs, vis_types, blend_types, save_float)
        # todo: add save_vrt=True to arguments and create a checkbox for this

    button_run_adaf.disabled = False


button_run_adaf.on_click(on_button_clicked)

# DISPLAY OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
output_widget = widgets.Output()

# DISPLAY WIDGET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~11
selection_options = widgets.HBox([
    widgets.VBox([
        widgets.Label("Select visualizations:"),
        dropdown_options_visualizations
    ]),
    widgets.VBox([
        widgets.Label("Select blends:"),
        dropdown_options_blends
    ])
])
display(
    widgets.VBox(
        [
            widgets.HTML(value=f"<b>RVT - TILED PROCESSING FOR LARGE FILES</b>"),
            b_file_select,
            selection_options,
            save_float_checkbox,
            button_run_adaf,
            output_widget
        ]
    )
)
