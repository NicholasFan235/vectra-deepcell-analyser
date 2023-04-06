# vectra-deepcell-analyser
Prototype for segmenting large immunofluorescence microscopy images using deepcell.
Deepcell works best with smaller images (up to 2000x2000 pixels).
This toolbox allows tiling large images, segmenting each tile with Deepcell and re-stitching the segmentations at the end.

# Summary:
Each component of the pipeline takes input files and creates output files for the next stages of the pipeline. (Good for prototyping but hard disk inefficient)


# Configuration
Tile size and deepcell parameters are in `config.deepcell_config.py`. DeepcellConfig class is used as a global singleton, so the parameters can be easily edited.

To take into account varying fluorescent stains, it is intended that a `Panel` class containing a summary of the stains is created.
Example see `panel_data/immune_panel.py`. This used for naming files and specifying which channels to give to deepcell.


# Running the pipeline
Input qptiffs should be placed in an input directory, structure like: `qptiffs/<EXPERIMENT_NAME>/<SAMPLE_NAME>.qptiff`
See example.py for running the code.


