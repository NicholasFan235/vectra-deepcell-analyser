import tifffile
import pathlib
import warnings

from ._helpers._get_files import _get_files
from .config import DeepcellConfig

def tile_for_deepcell(folder, name=None):
    if name is None:
        for file in _get_files(pathlib.Path('unstacked', folder), '.*\.tif'):
            tile_for_deepcell(folder, file.rstrip('.tif'))
    else:
        worker = _TileForDeepcell(folder, name)
        worker.process()


class _TileForDeepcell:
    def __init__(self, folder, name):
        self.folder = folder
        self.name = name

    def process(self):
        pathlib.Path('tiled_for_deepcell', self.folder, self.name).mkdir(exist_ok=True, parents=True)
        self._tile(tifffile.imread(
            pathlib.Path('unstacked', self.folder, f'{self.name}.tif')))
    
    def _tile(self, image):
        if image.ndim < 2:
            raise ValueError(f'Image should have ndim>=2, not {image.ndim}')

        for x in range(0, image.shape[1], DeepcellConfig.tile_width):
            for y in range(0, image.shape[0], DeepcellConfig.tile_height):
                x0 = max(0, x - DeepcellConfig.tile_padding_x)
                y0 = max(0, y - DeepcellConfig.tile_padding_y)
                x1 = min(x + DeepcellConfig.tile_width + DeepcellConfig.tile_padding_x, image.shape[1])
                y1 = min(y + DeepcellConfig.tile_height + DeepcellConfig.tile_padding_y, image.shape[0])
                tifffile.imwrite(
                    pathlib.Path('tiled_for_deepcell', self.folder, self.name,
                        f'{self.name}_{x0}_{y0}.png'),
                    image[y0:y1, x0:x1, ...])

