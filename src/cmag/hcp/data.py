# -*- coding: utf-8 -*-
###############################################################################
# hcp/data.py

"""Code and functions for loading and caching HCP subject data related to the
calculation of cortical magnification.

Because fitting models of cortical magnification to all 163 labeled HCP
subjects requires holding certain cortical-magnification-related information
for each subject in memory, we have to be conservative about using memory. We
also want to save time parsing the relevant parts of the data out of the HCP
data structures. For these reasons, the data loaded by this module is cached in
a directory. This allows one to later quickly load all the subject data without
parsing it out of the HCP data structures. The `cmag.hcp.config.cache_path`
value determines the default cache path for the data.

This module is intended to be used as a dataset object directly. The `load`
function can be used to load subject data.
"""


def load(sid, h,
         label_path=Ellipsis,
         cache_path=Ellipsis,
         overwrite=False):
    """Returns a dictionary of pRF and surface area data for the given HCP
    subject and hemisphere.

    This function loads HCP data using the `neuropythy` library, which must be
    properly configured (see https://github.com/noahbenson/neuropythy). Because
    loading data in this manner is expensive in terms of time and memory, the
    `load` function also generates cache files in a directory specified by the
    option `cache_path` and loads these files on subsequent invocations. The
    `cache_path` may be set to `Ellipsis` in order to use the cache path
    specified in the `cmag.hcp.config` subpackage.
    
    Parameters
    ----------
    sid : int
        The subject ID of the HCP subject whose data is to be loaded.
    h : 'lh' or 'rh'
        The hemisphere whose data is to be loaded.
    label_path : Ellipsis or path-like, optional
        The path from which the labels for the visual areas (i.e., V1, V2, V3,
        hV4, VO1, VO2 and V3a, V3b, IPS0, LO1) are found. The default value
        (`Ellipsis`) indicates that the path in `cmag.hcp.config.label_path`
        should be used.
    cache_path : Ellipsis or path-like, optional
        The directory into which cached files for the subject's data should be
        saved (or from which the data should be loaded). The default value
        (`Ellipsis`) indicates that the path in `cmag.hcp.config.cache_path`
        should be used. If `None` is given, then no caching is performed.
    overwrite : boolean, optional
        Whether to overwrite existing cache files, if found. The default is
        `False`.
    
    Returns
    -------
    dict
        A dictionary of data related to the subject/hemisphere. The dictionary
        contains the following keys: `'label', 'x', 'y', 'sigma',
        'polar_angle', 'eccentricity', 'surface_area', 'cod', 'midgray_x',
        'midgray_y', 'midgray_z', 'faces'`. All of these except for `faces` are
        vertex-wise values. The `faces` entry is a `3 x M` matrix of indices of
        the corner vertices for each of the triangle faces included in the
        visual areas.
    """
    import numpy as np
    import neuropythy as ny
    from pathlib import Path
    if cache_path is None:
        cachefile = None
    else:
        if cache_path is Ellipsis:
            from .config import cache_path
        vfile = Path(cache_path) / f'{h}.{sid}_vert.npy'
        ffile = Path(cache_path) / f'{h}.{sid}_face.npy'
    if overwrite or not (cache_path and vfile.is_file() and ffile.is_file()):
        if label_path is Ellipsis:
            from .config import label_path
        sub = ny.data['hcp_lines'].subjects[sid]
        hem = sub.hemis[h]
        lbl = np.array(hem.prop('visual_area'), dtype=np.int16)
        labelfile = Path(label_path.format(sid=sid, hemisphere=h))
        if not labelfile.is_file():
            raise RuntimeError(f"label file not found for subject {sid}/{h}")
        lbl = ny.load(str(labelfile))
        nz = (lbl > 0)
        msh = hem.surface('midgray').submesh(nz)
        vid = msh.labels
        prf_x = msh.prop('prf_x')
        prf_y = msh.prop('prf_y')
        prf_sig = msh.prop('prf_radius')
        prf_ang = msh.prop('prf_polar_angle')
        prf_ecc = msh.prop('prf_eccentricity')
        prf_cod = msh.prop('prf_variance_explained')
        surfarea = msh.prop('midgray_surface_area')
        mid_x = msh.coordinates[0]
        mid_y = msh.coordinates[1]
        mid_z = msh.coordinates[2]
        vdata = np.stack(
            [vid, lbl[vid],
             prf_x, prf_y, prf_sig, prf_ang, prf_ecc, prf_cod,
             surfarea, mid_x, mid_y, mid_z],
            axis=0,
            dtype=np.float32)
        faces = msh.tess.indexed_faces
        if vfile:
            np.save(vfile, vdata)
        if ffile:
            np.save(ffile, faces)
    else:
        vdata = np.load(vfile)
        faces = np.load(ffile)
    vdata = vdata.astype(float)
    return dict(
        vertex_label=vdata[0].astype(np.int32),
        label=vdata[1].astype(np.int32),
        x=vdata[2], y=vdata[3], sigma=vdata[4],
        polar_angle=vdata[5], eccentricity=vdata[6], cod=vdata[7],
        surface_area=vdata[8],
        coordinates=vdata[9:],
        faces=faces)

def joinhemis(lhdata, rhdata=None, /):
    """Joins left and right hemisphere HCP data into a single data object.

    The `cmag.hcp.data.load` function returns a data structure (a dictionary of
    numpy arrays) for an individual hemisphere. To run analyses on the entire
    bilateral visual area, these hemispheres must be joined. This function can
    be called either as `joinhemis(lhdata, rhdata)` or as 
    `joinhemis((lhdata, rhdata))`. It returns the data with all keys joined and
    the additional key `'hemisphere'` whose values are either `'lh'` or `'rh'`.
    """
    from numpy import full, concatenate
    if rhdata is None:
        (lhdata, rhdata) = lhdata
    basic_keys = (
        'vertex_label', 'label', 'x', 'y', 'sigma', 'polar_angle',
        'eccentricity', 'cod', 'surface_area', 'coordinates', 'faces')
    data = {
        k: concatenate([lhdata[k], rhdata[k]], axis=-1)
        for k in basic_keys}
    nlh = len(lhdata['label'])
    nrh = len(rhdata['label'])
    data['hemisphere'] = concatenate([full(nlh, 'lh'), full(nrh, 'rh')])
    return data
