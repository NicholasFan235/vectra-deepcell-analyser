import tifffile
import xmltodict
import pathlib
import warnings

from ._helpers._get_files import _get_files
from .panel_data import ImmunePanel


def split_immune_qptiff(folder, name=None):
    if name is None:
        for file in _get_files(pathlib.Path('qptiffs', folder), '.*\.qptiff'):
            split_immune_qptiff(folder, file.rstrip('.qptiff'))
    else:
        worker = _SplitQpTiff(folder, name, ImmunePanel)
        worker.process()


class _SplitQpTiff:
    def __init__(self, folder, name, config):
        self.folder = folder
        self.name = name
        self.config = config

    def process(self):
        input_file = pathlib.Path('qptiffs', self.folder, f'{self.name}.qptiff')
        if not input_file.exists() or not input_file.is_file():
            raise FileNotFoundError(f'{input_file} was not found or is a directory')
        tf = tifffile.TiffFile(input_file)
        pathlib.Path('unstacked', self.folder).mkdir(exist_ok=True, parents=True)
        self._process_tiff(tf)

    def _process_tiff(self, tf:tifffile.TiffFile):
        for i in range(self.config.n_channels):
            self._process_page(tf.pages[i])

    def _process_page(self, page):
        assert page.tags[270].name == "ImageDescription",\
            f"Expected tag 270: ImageDescription, found {page.tags[270].name}"
        im_desc = xmltodict.parse(page.tags[270].value)
        channel = im_desc['PerkinElmer-QPI-ImageDescription']['Name']
        if not channel in self.config.channel_map:
            warnings.warn(f"Unmapped channel \"{channel}\", ignoring.", UserWarning)
            return
        
        tifffile.imwrite(
            pathlib.Path('unstacked', self.folder, f"{self.name}_{self.config.channel_map[channel]}.tif"),
            page.asarray())
