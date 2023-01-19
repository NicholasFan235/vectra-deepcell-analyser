import tifffile
import numpy as np
import pathlib
import tqdm
import pandas as pd
import geojson
from .config import DeepcellConfig

def calculate_centroids(folder, name,
        compartment='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    worker = _CalculateCetroidsWorker(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    worker.process()

def centroid_geojson(folder, name,
        compartment='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    pass


class _CalculateCetroidsWorker:
    def __init__(self, folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold):
        self.folder = folder
        self.name = name
        self.compartment = compartment
        self.nucleus_channel = nucleus_channel
        self.membrane_channel = membrane_channel
        self.interior_threshold = interior_threshold if interior_threshold is not None else DeepcellConfig.interior_threshold
        self.maxima_threshold = maxima_threshold if maxima_threshold is not None else DeepcellConfig.maxima_threshold

        self.infile_basename = f'{self.name}_{self.compartment}_{self.nucleus_channel}_{self.membrane_channel}_{str(int(self.interior_threshold*1000))}_{str(int(self.maxima_threshold*1000))}'

    def has_processed(self):
        outfolder = pathlib.Path('centroids', self.folder, self.name)
        outfolder.mkdir(exist_ok=True, parents=True)
        return pathlib.Path(outfolder,  f'{self.infile_basename}_centroids.csv').is_file()

    def process(self):
        infile = pathlib.Path('deepcell_labelled', self.folder, self.name, f'{self.infile_basename}.tif')
        if not infile.is_file():
            raise FileNotFoundError(f'{infile} does not exist or is a directory')
        
        outfolder = pathlib.Path('centroids', self.folder, self.name)
        outfolder.mkdir(exist_ok=True, parents=True)

        data = pd.DataFrame(columns=('Object Id', 'centroid_x_pixels', 'centroid_y_pixels')).set_index('Object Id')

        im = tifffile.imread(infile)
        labelled = im.reshape(im.shape + (1,))
        indices = np.moveaxis(np.indices(im.shape), 0, 2)
        formatted = np.dstack((indices, labelled))[im>0,...].reshape((-1, 3))

        df = pd.DataFrame(formatted, columns=['y', 'x', 'label'])
        df.to_csv(pathlib.Path(outfolder, f'{self.infile_basename}_pixel_data.csv'), index=False)

        centroids = df.groupby('label').mean()
        centroids = centroids.rename(columns={'y':'centroid_y_pixels', 'x':'centroid_x_pixels'})
        centroids['area_pixels'] = df.groupby('label').size()
        centroids = centroids.join(df.groupby('label').min().rename(columns={'y':'min_y_pixels', 'x':'min_x_pixels'}))
        centroids = centroids.join(df.groupby('label').max().rename(columns={'y':'max_y_pixels', 'x':'max_x_pixels'}))

        centroids.index = centroids.index.rename('Object Id')
        centroids.to_csv(pathlib.Path(outfolder, f'{self.infile_basename}_centroids.csv'))

    def make_geojson(self):
        outfolder = pathlib.Path('centroids', self.folder, self.name)
        outfolder.mkdir(exist_ok=True, parents=True)
        if not pathlib.Path(outfolder,  f'{self.infile_basename}_centroids.csv').is_file():
            self.process()
        outfile = pathlib.Path(outfolder, f'{self.infile_basename}_centroids.geojson')

        centroids = pd.read_csv(pathlib.Path(outfolder, f'{self.infile_basename}_centroids.csv'), index_col='Object Id')
        points = geojson.MultiPoint(
            list(centroids[['centroid_x_pixels', 'centroid_y_pixels']].itertuples(index=False, name=None)),
            properties={"object_type": "detection", "isLocked": True, "Name": "deepcell"})
        
        geojson.dump([points], open(outfile, 'w'), indent=2)
