import tifffile
import pathlib
import pandas as pd
import numpy as np

from .config import DeepcellConfig
from .panel_data import ImmunePanel


def compute_immune_markers(folder, name,
        compartment:str='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    worker = _CalculateMarkers(folder, name, ImmunePanel, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    worker.process()


class _CalculateMarkers:
    def __init__(self, folder, name, panel, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold):
        self.folder = folder
        self.name = name
        self.panel = panel
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
        self.unstacked_folder = pathlib.Path('unstacked', self.folder)
        for marker in self.panel.channel_map.values():
            p = pathlib.Path(self.unstacked_folder, f'{self.name}_{marker}.tif')
            if not p.is_file():
                raise FileNotFoundError(f'{p} does not exist or is a directory')

    def process(self):
        arrs = []
        markers = ['y', 'x', 'label']
        im = tifffile.imread(self.labelled_file)
        arrs.append(im.reshape(im.shape + (1,)))
        arrs.append(np.moveaxis(np.indices(im.shape), 0, 2))

        for marker in self.panel.channel_map.values():
            markers.append(marker)
            arrs.append(self._get_marker(marker))

        formatted = np.dstack(arrs)[im>0, ...].reshape((-1, len(markers)))
        return pd.DataFrame(formatted, columns=markers).set_index(['y', 'x'])

    def _get_marker(self, marker):
        im = tifffile.imread(pathlib.Path(self.unstacked_folder, f'{self.name}_{marker}.tif'))
        return im.reshape(im.shape + (1,))
