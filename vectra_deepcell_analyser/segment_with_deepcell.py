import tifffile
import pathlib
import re
import tqdm
import os
import numpy as np

from .config import DeepcellConfig

def segment_with_deepcell(
        folder, name,
        compartment='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    worker = _DeepcellWorker(
        folder, name,
        compartment, nucleus_channel, membrane_channel,
        interior_threshold, maxima_threshold)
    worker.process()


class _DeepcellWorker:
    def __init__(self, folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold):
        self.folder = folder
        self.name = name
        self.compartment = compartment
        self.nucleus_channel = nucleus_channel
        self.membrane_channel = membrane_channel
        self.interior_threshold = interior_threshold if interior_threshold is not None else DeepcellConfig.interior_threshold
        self.maxima_threshold = maxima_threshold if maxima_threshold is not None else DeepcellConfig.maxima_threshold

        self.nuc_folder = pathlib.Path('tiled_for_deepcell', self.folder, f'{self.name}_{self.nucleus_channel}')
        if not self.nuc_folder.is_dir():
            raise FileNotFoundError(f'{self.nuc_folder} does not exist or is not a directory')
        self.mem_folder = pathlib.Path('tiled_for_deepcell', self.folder, f'{self.name}_{self.membrane_channel}')
        if not self.mem_folder.is_dir():
            raise FileNotFoundError(f'{self.mem_folder} does not exist or is not a directory')

    def process(self):
        import deepcell
        pathlib.Path('deepcell_labelled_tiles', self.folder, self.name).mkdir(
            exist_ok=True, parents=True)

        pattern = re.compile(f'{self.name}_{self.nucleus_channel}_(?P<x0>\d+)_(?P<y0>\d+)\.png')
        for file in os.listdir(self.nuc_folder):
            m = pattern.match(file)
            self._process_tile(int(m['x0']), int(m['y0']))

    def _process_tile(self, x0:int, y0:int):
        print(x0, y0)
        nuc = tifffile.imread(
            pathlib.Path(self.nuc_folder, f'{self.name}_{self.nucleus_channel}_{x0}_{y0}.png'))
        mem = tifffile.imread(
            pathlib.Path(self.mem_folder, f'{self.name}_{self.membrane_channel}_{x0}_{y0}.png'))
        im = np.stack((nuc/np.max(nuc), mem/np.max(mem)), axis=-1)
        im = np.expand_dims(im, 0)

        app = deepcell.applications.Mesmer()
        interior_threshold = self.interior_threshold
        maxima_threshold = self.maxima_threshold
        predictions = app.predict(im,
            image_mpp=DeepcellConfig.image_mpp,
            postprocess_kwargs_whole_cell={
                'interior_threshold': interior_threshold,
                'maxima_threshold': maxima_threshold},
            compartment = self.compartment)
        interior_threshold_str = str(int(interior_threshold*1000))
        maxima_threshold_str = str(int(maxima_threshold*1000))

        tifffile.imwrite(
            pathlib.Path('deepcell_labelled_tiles', self.folder, self.name,
                f'{self.name}_{self.compartment}_{self.nucleus_channel}_{self.membrane_channel}_{interior_threshold_str}_{maxima_threshold_str}',
                f'{self.name}_{self.compartment}_{self.nucleus_channel}_{self.membrane_channel}_{interior_threshold_str}_{maxima_threshold_str}_{x0}_{y0}.tif'),
            predictions[0,...])
