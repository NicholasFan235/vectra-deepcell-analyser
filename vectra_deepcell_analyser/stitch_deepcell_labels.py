import numpy as np
import pathlib
import os
import re
import tifffile

from .config import DeepcellConfig

def stitch_deepcell_labels(
        folder, name,
        compartment:str='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    x_stitch_worker = _StitchDeepcellLabelsX(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    x_stitch_worker.process()
    
    y_stitch_worker = _StitchDeepcellLabelsY(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    y_stitch_worker.process()

def stitch_deepcell_labels_y(
        folder, name,
        compartment:str='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    y_stitch_worker = _StitchDeepcellLabelsY(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    y_stitch_worker.process()

def stitch_deepcell_labels_x(
        folder, name,
        compartment:str='whole-cell',
        nucleus_channel:str='DAPI',
        membrane_channel:str='ECad',
        interior_threshold:float=None,
        maxima_threshold:float=None):
    x_stitch_worker = _StitchDeepcellLabelsX(folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold)
    x_stitch_worker.process()


class _StitchDeepcellLabelsX:
    def __init__(self, folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold):
        self.folder = folder
        self.name = name
        self.compartment = compartment
        self.nucleus_channel = nucleus_channel
        self.membrane_channel = membrane_channel
        self.interior_threshold = interior_threshold if interior_threshold is not None else DeepcellConfig.interior_threshold
        self.maxima_threshold = maxima_threshold if maxima_threshold is not None else DeepcellConfig.maxima_threshold
    
        p = pathlib.Path('unstacked', folder, f'{name}_{nucleus_channel}.tif')
        if not p.is_file():
            raise FileNotFoundError(f'{p} does not exist or is a directory. Used to check original image dimensions')
        im = tifffile.TiffFile(p)
        self.original_shape = im.pages[0].shape

        self.tile_basename = f'{name}_{compartment}_{nucleus_channel}_{membrane_channel}_{str(int(self.interior_threshold*1000))}_{str(int(self.maxima_threshold*1000))}'
        self.tiles_folder = pathlib.Path('deepcell_labelled_tiles', folder, name, self.tile_basename)
        if not self.tiles_folder.is_dir():
            raise FileNotFoundError(f'{self.tiles_folder} is not a directory')

    def process(self):
        self.outfolder = pathlib.Path('deepcell_labelled_tiles_x_stitched', self.folder, self.name, self.tile_basename)
        self.outfolder.mkdir(exist_ok=True, parents=True)

        for y in range(0, self.original_shape[0], DeepcellConfig.tile_height):
            y0 = max(0, y-DeepcellConfig.tile_padding_y)
            self._stitch_x(y0)

    def _stitch_x(self, y):
        print(f'Stitching y0:{y}')
        ims = []
        for x in range(0, self.original_shape[1], DeepcellConfig.tile_width):
            x0 = max(0, x - DeepcellConfig.tile_padding_x)
            im = tifffile.TiffFile(pathlib.Path(
                self.tiles_folder, f'{self.tile_basename}_{x0}_{y}.tif'))
            ims.append(im)
        stitched = np.zeros((ims[0].pages[0].shape[0], self.original_shape[1]), dtype='uint32')

        offsets = []
        for i, im in enumerate(ims):
            x0 = i*DeepcellConfig.tile_width
            x1 = (i+1)*DeepcellConfig.tile_width
            x0 += DeepcellConfig.tile_padding_x if i > 0 else 0
            x1 -= DeepcellConfig.tile_padding_x if i < len(ims)-1 else 0
            tile_x0 = 2*DeepcellConfig.tile_padding_x
            tile_x1 = tile_x0 + DeepcellConfig.tile_width - 2*DeepcellConfig.tile_padding_x
            tile_x0 = 0 if i <= 0 else tile_x0
            tile_x1 -= DeepcellConfig.tile_padding_x if i <= 0 else 0
            
            offsets.append(np.max(stitched))
            tile = im.pages[0].asarray()[:, tile_x0:tile_x1]
            tile[tile>0] += offsets[i]
            stitched[:, x0:x1] = tile
        print(f'Done basic stitch y0:{y}')

        for i in range(len(ims)-1):
            print(f'Solving overlap #{i}')
            x0 = (i+1)*DeepcellConfig.tile_width - DeepcellConfig.tile_padding_x
            x1 = (i+1)*DeepcellConfig.tile_width + DeepcellConfig.tile_padding_x
            lx0 = DeepcellConfig.tile_width
            if i == 0: lx0 -= DeepcellConfig.tile_padding_x
            rx1 = DeepcellConfig.tile_padding_x * 2
            l = ims[i].pages[0].asarray()[:, lx0:]
            l[l>0] += offsets[i]
            r = ims[i+1].pages[0].asarray()[:, :rx1]
            r[r>0] += offsets[i+1]
            self._solve_overlap(
                ims[i].pages[0].asarray()[:,lx0:],
                ims[i+1].pages[0].asarray()[:,:rx1],
                stitched, x0, x1)
        
        outfile = pathlib.Path(self.outfolder, f'{self.tile_basename}_{y}.tif')
        tifffile.imwrite(outfile, stitched)

    def _solve_overlap(self, l, r, stitched, x0, x1):
        cs = set(l[:,l.shape[1]//2].reshape(-1))
        cs.add(0)
        cs.remove(0)
        lmask = np.zeros(l.shape, dtype=bool)
        lmask[:,:l.shape[1]//2] = True
        overlap_mask = np.isin(l, list(cs))
        lmask[overlap_mask] = True

        rmask = np.ones(l.shape, dtype=bool)
        rmask[:,:l.shape[1]//2] = False
        cs = set(r[:,r.shape[1]//2]).union(set(r[overlap_mask]))
        cs.add(0)
        cs.remove(0)
        rmask[np.isin(r, list(cs))] = False

        stitched[:,x0:x1][lmask] = l[lmask]
        stitched[:,x0:x1][rmask] = r[rmask]
        stitched = self._fill_holes(l, r, lmask, rmask, stitched, x0, x1)
        return stitched

    def _fill_holes(self, l, r, lmask, rmask, stitched, x0, x1):
        print(f'Filling holes')
        candidates = (~(lmask|rmask))&((l>0)|(r>0))
        missing = 0
        for i, j in zip(*candidates.nonzero()):
            if stitched[:,x0:x1][i, j] > 0 or l[i,j]<=0 or r[i,j]<=0:
                continue
            lm = (l==l[i,j]) & (stitched[:,x0:x1]<=0)
            rm = (r==r[i,j]) & (stitched[:,x0:x1]<=0)
            m = lm|rm
            if m.sum() > 50:
                missing += 1
                stitched[:,x0:x1][m] = np.max(stitched)+1
        print(f'Filled {missing} missing cells')
        return stitched


class _StitchDeepcellLabelsY:
    def __init__(self, folder, name, compartment, nucleus_channel, membrane_channel, interior_threshold, maxima_threshold):
        self.folder = folder
        self.name = name
        self.compartment = compartment
        self.nucleus_channel = nucleus_channel
        self.membrane_channel = membrane_channel
        self.interior_threshold = interior_threshold if interior_threshold is not None else DeepcellConfig.interior_threshold
        self.maxima_threshold = maxima_threshold if maxima_threshold is not None else DeepcellConfig.maxima_threshold
    
        p = pathlib.Path('unstacked', folder, f'{name}_{nucleus_channel}.tif')
        if not p.is_file():
            raise FileNotFoundError(f'{p} does not exist or is a directory. Used to check original image dimensions')
        im = tifffile.TiffFile(p)
        self.original_shape = im.pages[0].shape

        self.tile_basename = f'{name}_{compartment}_{nucleus_channel}_{membrane_channel}_{str(int(self.interior_threshold*1000))}_{str(int(self.maxima_threshold*1000))}'
        self.tiles_folder = pathlib.Path('deepcell_labelled_tiles_x_stitched', folder, name, self.tile_basename)
        if not self.tiles_folder.is_dir():
            raise FileNotFoundError(f'{self.tiles_folder} is not a directory')

    def process(self):
        self.outfolder = pathlib.Path('deepcell_labelled', self.folder, self.name)
        self.outfolder.mkdir(exist_ok=True, parents=True)

        self._stitch_y()

    def _stitch_y(self):
        ims = []
        for y in range(0, self.original_shape[0], DeepcellConfig.tile_height):
            y0 = max(0, y - DeepcellConfig.tile_padding_y)
            im = tifffile.TiffFile(pathlib.Path(
                self.tiles_folder, f'{self.tile_basename}_{y0}.tif'))
            ims.append(im)
        stitched = np.zeros(self.original_shape, dtype='uint32')

        offsets = []
        for i, im in enumerate(ims):
            y0 = i*DeepcellConfig.tile_height
            y1 = (i+1)*DeepcellConfig.tile_height
            y0 += DeepcellConfig.tile_padding_y if i > 0 else 0
            y1 -= DeepcellConfig.tile_padding_y if i < len(ims)-1 else 0
            tile_y0 = 2*DeepcellConfig.tile_padding_y
            tile_y1 = tile_y0 + DeepcellConfig.tile_height - 2*DeepcellConfig.tile_padding_y
            tile_y0 = 0 if i <= 0 else tile_y0
            tile_y1 -= DeepcellConfig.tile_padding_y if i <= 0 else 0
            
            offsets.append(np.max(stitched))
            tile = im.pages[0].asarray()[tile_y0:tile_y1,:]
            tile[tile>0] += offsets[i]
            stitched[y0:y1,:] = tile

        for i in range(len(ims)-1):
            outfile = pathlib.Path(self.outfolder, f'{self.tile_basename}.tif')
            tifffile.imwrite(outfile, stitched)
            print(f'Solving overlap #{i}')
            y0 = (i+1)*DeepcellConfig.tile_height - DeepcellConfig.tile_padding_y
            y1 = (i+1)*DeepcellConfig.tile_height + DeepcellConfig.tile_padding_y
            ly0 = DeepcellConfig.tile_height
            if i == 0: ly0 -= DeepcellConfig.tile_padding_y
            ry1 = DeepcellConfig.tile_padding_y * 2
            l = ims[i].pages[0].asarray()[ly0:, :]
            l[l>0] += offsets[i]
            r = ims[i+1].pages[0].asarray()[:ry1, :]
            r[r>0] += offsets[i+1]
            self._solve_overlap(
                l, r, stitched, y0, y1)
        
        outfile = pathlib.Path(self.outfolder, f'{self.tile_basename}.tif')
        tifffile.imwrite(outfile, stitched)

    def _solve_overlap(self, l, r, stitched, y0, y1):
        cs = set(l[l.shape[0]//2,:].reshape(-1))
        cs.add(0)
        cs.remove(0)
        lmask = np.zeros(l.shape, dtype=bool)
        lmask[:l.shape[0]//2,:] = True
        overlap_mask = np.isin(l, list(cs))
        lmask[overlap_mask] = True

        rmask = np.ones(l.shape, dtype=bool)
        rmask[:l.shape[0]//2,:] = False
        cs = set(r[r.shape[0]//2,:]).union(set(r[overlap_mask]))
        cs.add(0)
        cs.remove(0)
        rmask[np.isin(r, list(cs))] = False

        stitched[y0:y1,:][lmask] = l[lmask]
        stitched[y0:y1,:][rmask] = r[rmask]
        stitched = self._fill_holes(l, r, lmask, rmask, stitched, y0, y1)
        return stitched

    def _fill_holes(self, l, r, lmask, rmask, stitched, y0, y1):
        print(f'Filling holes')
        candidates = (~(lmask|rmask))&((l>0)|(r>0))
        missing = 0
        for i, j in zip(*candidates.nonzero()):
            if stitched[y0:y1,:][i, j] > 0 or l[i,j]<=0 or r[i,j]<=0:
                continue
            lm = (l==l[i,j]) & (stitched[y0:y1,:]<=0)
            rm = (r==r[i,j]) & (stitched[y0:y1,:]<=0)
            m = lm|rm
            if m.sum() > 50:
                missing += 1
                stitched[y0:y1,:][m] = np.max(stitched)+1
        print(f'Filled {missing} missing cells')
        return stitched