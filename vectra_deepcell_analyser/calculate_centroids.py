import tifffile
import numpy as np
import pathlib
import tqdm
import pandas as pd
from .config import DeepcellConfig

def calculate_centroids(folder, name,
        compartment='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    worker = _CalculateCetroidsWorker(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    worker.process()


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


    def process(self):
        infile = pathlib.Path('deepcell_labelled', self.folder, f'{self.infile_basename}.tif')
        if not infile.is_file():
            raise FileNotFoundError(f'{infile} does not exist or is a directory')
        
        outfolder = pathlib.Path('centroids', self.folder)
        outfolder.mkdir(exist_ok=True, parents=True)

        data = pd.DataFrame(columns=('Object Id', 'centroid_x_pixels', 'centroid_y_pixels')).set_index('Object Id')

        labelled = tifffile.imread(infile)
        max_cid = np.max(labelled)
        for id in tqdm.tqdm(range(max_cid)):
            m = np.argwhere(labelled==id)
            if m.shape[0] > 0:
                data.loc[id] = m.mean(0)

        data['centroid_x'] = data.centroid_x_pixels * DeepcellConfig.image_mpp
        data['centroid_y'] = data.centroid_y_pixels * DeepcellConfig.image_mpp

        data.to_csv(pathlib.Path('centroids', self.folder, f'{self.infile_basename}_centroids.csv'))
