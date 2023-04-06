import vectra_deepcell_analyser as vda
# timer for logging
from vectra_deepcell_analyser._helpers._timer import _Timer


# Configure
vda.config.DeepcellConfig.tile_height = 4000
vda.config.DeepcellConfig.tile_width = 4000
vda.config.DeepcellConfig.tile_padding_x = 100
vda.config.DeepcellConfig.tile_padding_y = 100

folder = 'test'
# If file is None, will run on all files in folder
file = 'TEST_SAMPLE'


# Run
'''
Split the original qptiff into multiple .tif images, one for each channel.
Files are named based on the Panel class given.
This step expects qptiffs from vectra and so expects a piece of metadata describing the wavelength in the qptiff stack
Currently specific ImmunePanel is hard coded in, TODO: create more general method

Input: 'qptiffs/<folder>/<file>.qptiff'
Outputs: 'unstacked/<folder>/<file>_<marker>.tif'
'''
with _Timer('Split'):
    vda.split_immune_qptiff(folder, file)

'''
Tiles the split images based on configured DeepcellConfig static variables

Input: 'unstacked/<folder>/<file>_*.tif
Output: 'tiled_for_deepcell/<folder>/<file>_<marker>/<file>_<marker>_<x>_<y>.png'
'''
with _Timer('Tile'):
    vda.tile_for_deepcell(folder, file)

'''
Averages all markers which are not DAPI
This can be useful in case there is no specific membrane marker.
In this case, the average of all other markers is used to approximate a membrane marker. 
This works on the tiled images to avoid memory issues

Input: 'tiled_for_deepcell/<folder>/<file>_<marker>/<file>_<marker>_<x>_<y>.png'
Output: 'tiled_for_deepcell/<folder>/<file>_AVGMARKER/<file>_AVGMARKER_<x>_<y>.png'
'''
with _Timer('Mean Membrane Marker'):
    vda.generate_mean_immune_marker(folder, file)

'''
Segment each tile using deepcell.
Specify Deepcell compartment here, defaults to 'whole-cell'
Deepcell parameters can be passed in as arguments here, otherwise parameters set in DeepcellConfig will be used.
Also specify which channels to use as nuclear and membrane markers. Defaults to find channels called DAPI and ECad.

Input: 'tiled_for_deepcell/<folder>/<file>_<membrane marker>/<file>_<membrane marker>_<x>_<y>.png'
       'tiled_for_deepcell/<folder>/<file>_<nuclear marker>/<file>_<nuclear marker>_<x>_<y>.png'
Output: 'deepcell_labelled_tiles/<folder>/<file>/<file>_<deepcell config>/<file>_<deepcell config>_<x>_<y>.tif'
'''
with _Timer('Segment'):
    vda.segment_with_deepcell(folder, file)

'''
These tiles are then stitched back together.
Note currently this is very slow as stitch goes through many steps.
Tiles are first stitched along x-direction and written to file.
Once everyting is stitched along x-axis, stitch together along y.

Input: 'deepcell_labelled_tiles/<folder>/<file>/<file>_<deepcell config>/<file>_<deepcell config>_<x>_<y>.tif'
Intermediate: 'deepcell_labelled_tiles_x_stitched/<folder>/<file>/<file>_<deepcell config>/<file>_<deepcell config>_<y>.tif'
Output: 'deepcell_labelled/<folder>/<file>/<file>_<deepcell config>.tif'
'''
with _Timer('Stitch'):
    vda.stitch_deepcell_labels(folder, file)


# Postprocessing
'''
Once segmented may want to convert segmentation to cell outlines
Currently this is quite clunky, since can't add as overlay to original image.
Instead create new image from nuclear and membrane markers, then draw on segmented cell outlines.
Therefore have to provide nuclear and membrane channels, as well as used deepcell config to this method

Input: 'deepcell_labelled/<folder>/<file>/<file>_<deepcell config>.tif'
       'unstacked/<folder>/<file>/<file>_<membrane marker>.tif'
       'unstacked/<folder>/<file>/<file>_<nuclear marker>.tif'
Output: 'outlined/<folder>/<file>_<outlined>.tif'
'''
# Can easily create memory issues since loading huge images to RAM
#with _Timer('Make Outline Overlay'):
#    vda.make_outline_overlay(folder, file)

'''
Useful to get centroids from segmentation
Note currently very memory inefficient, but not had problems running on clusters (probably won't run on laptop)
Creates csv with cell centroids
Also creates geojson which can be loaded into visualisation software. Tested with QuPath, to overlay centroids onto original image

Input: 'deepcell_labelled/<folder>/<file>/<file>_<deepcell config>.tif'
Output: 'centroids/<folder>/<file>/<file>_<deepcell config>.csv'
        'centroids/<folder>/<file>/<file>_<deepcell config>.geojson'
'''
with _Timer('Calculate Centroids'):
    vda.calculate_centroids(folder, file)

'''
Can also use segmentation and calculate average marker intensity across each cell.
Requires Panel to be specified (so it knows what channels to look for)
Currently only for specific ImmunePanel, TODO: create general version of method for new panels
Creates csv with cell id, centroid and average marker intensities

Input: 'deepcell_labelled/<folder>/<file>/<file>_<deepcell config>.tif'
       'unstacked/<folder>/<file>/<file>_<markers>.tif'
Output: 'output/<folder>/<file>/<file>_<deepcell config>.csv'
'''
with _Timer('Calculate Avg Marker Intenstity For Each Cell'):
    vda.compute_immune_markers(folder, file)

    