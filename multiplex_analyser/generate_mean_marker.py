import tifffile
import numpy as np
import pathlib
import os
import re
import tqdm

from .panel_data import ImmunePanel
from .config import DeepcellConfig


def generate_mean_immune_marker(folder, name):
    worker = _GenerateMeanMarker(folder, name, ImmunePanel)
    worker.process()


class _GenerateMeanMarker:
    def __init__(self, folder, name, panel):
        self.folder = folder
        self.name = name
        self.panel = panel

    def process(self):
        dapi_tile_folder = pathlib.Path('tiled_for_deepcell', self.folder, f'{self.name}_DAPI')
        if not dapi_tile_folder.is_dir():
            raise FileNotFoundError('{dapi_tile_folder} does not exist or is not a directory')

        pathlib.Path('tiled_for_deepcell', self.folder, f'{self.name}_AVGMARKER').mkdir(
            exist_ok=True, parents=True)
        
        pattern = re.compile(f'{self.name}_DAPI_(?P<x0>\d+)_(?P<y0>\d+)\.png')
        for file in tqdm.tqdm(os.listdir(dapi_tile_folder)):
            m = pattern.match(file)
            self._process_tile(int(m['x0']), int(m['y0']))

        
    def _process_tile(self, x0:int, y0:int):
        mean_im = None
        n_channels = 0
        
        for _, channel in self.panel.channel_map.items():
            if channel == "DAPI":
                continue
            new_im = self._get_channel(x0, y0, channel)

            if mean_im is None:
                mean_im = new_im
            else:
                mean_im = mean_im + new_im
            n_channels += 1
        
        mean_im = mean_im / n_channels
        mean_im = mean_im.astype('uint8')
        tifffile.imwrite(
            pathlib.Path('tiled_for_deepcell', self.folder,
                f'{self.name}_AVGMARKER',
                f'{self.name}_AVGMARKER_{x0}_{y0}.png'),
            mean_im)

    def _get_channel(self, x0:int, y0:int, channel:str):
        image_file = pathlib.Path('tiled_for_deepcell', self.folder,
            f'{self.name}_{channel}',
            f'{self.name}_{channel}_{x0}_{y0}.png')
        if not image_file.is_file() or not image_file.exists():
            raise FileNotFoundError(f'{image_file} was not found or is a directory')
        im = tifffile.imread(image_file)
        if not im.ndim == 2:
            raise ValueError(f'{image_file} should have ndim=2, not {im.ndim}')
        return im.astype(int)
