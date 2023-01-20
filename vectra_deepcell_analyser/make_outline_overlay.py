import tifffile
import skimage.segmentation
import pathlib
import numpy as np

from .config import DeepcellConfig


def make_outline_overlay(folder, name,
        compartment:str='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    worker = _OutlineOverlayWorker(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    worker.process()
    

class _OutlineOverlayWorker:
    def __init__(self, folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold):
        self.folder = folder
        self.name = name
        self.compartment = compartment
        self.nucleus_channel = nucleus_channel
        self.membrane_channel = membrane_channel
        self.interior_threshold = interior_threshold if interior_threshold is not None else DeepcellConfig.interior_threshold
        self.maxima_threshold = maxima_threshold if maxima_threshold is not None else DeepcellConfig.maxima_threshold

        self.deepcell_basename = f'{name}_{compartment}_{nucleus_channel}_{membrane_channel}_{str(int(self.interior_threshold*1000))}_{str(int(self.maxima_threshold*1000))}'
        self.labelled_file = pathlib.Path(
            'deepcell_labelled', self.folder, self.name, f'{self.deepcell_basename}.tif')
        if not self.labelled_file.is_file():
            raise FileNotFoundError(f'{self.labelled_file} does not exist or is a directory')
        self.nucleus_file = pathlib.Path(
            'unstacked', self.folder, f'{self.name}_{self.nucleus_channel}.tif')
        if not self.nucleus_file.is_file():
            raise FileNotFoundError(f'{self.nucleus_file} does not exist or is a directory')
        self.membrane_file = pathlib.Path(
            'unstacked', self.folder, f'{self.name}_{self.membrane_channel}.tif')
        if not self.membrane_file.is_file():
            raise FileNotFoundError(f'{self.membrane_file} does not exist or is a directory')

    def process(self):
        outfolder = pathlib.Path('outlined', self.folder)
        outfolder.mkdir(exist_ok=True, parents=True)
        outfile = pathlib.Path(self.outfolder, f'{self.deepcell_basename}_outlined.qptiff')

        labelled = tifffile.imread(self.labelled_file)
        outlines = skimage.segmentation.find_boundaries(labelled, connectivity=1, mode='inner')
        
        nuc = tifffile.imread(self.nucleus_file)
        mem = tifffile.imread(self.membrane_file)

        rgb_image = np.stack((mem*0, mem/np.max(mem),nuc/np.max(nuc)), axis=-1)
        rgb_image = np.expand_dims(rgb_image, 0)

        rgb_image[outlines>0, :] = 1
        tifffile.imwrite(outfile, rgb_image)
