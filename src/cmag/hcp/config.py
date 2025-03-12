# -*- coding: utf-8 -*-
###############################################################################
# hcp/config.py

"""Configuration variables for the HCP analysis component of the cmag library.

The attributes of this subpackage of `cmag.hcp` are intended to be changed upon
first loading the library in order to match one's system configuration.
Specifically, the `label_path` and `cache_path` attributes should be updated.
Other functions in the `hcp` subpackage that accept an option named after one
of the attributes in this subpackage (e.g., `label_path` or `cache_path`) will
convert an argument of `Ellipsis` (`...`) into the current value from this
subpackage. (In other words, if function `f` takes an option `cache_path` and
it is invoked as `f(*args, cache_path=Ellipsis)`, then `f` will import the
current value of `cache_path` from `cmag.hcp.config` and use it.) Note that
most functions use `Ellipsis` as the default value for these options.

The `label_names` attribute should not be updated.
"""


# Utilities ###################################################################

def getenv(name, default=None):
    """Returns the value associated with the given OS environment variable.

    If the name is not found then the optional parameter `default` is returned.
    If not provided, this value is `None`. Note that this function never raises
    an error in the case that a value is not found.
    """
    from os import environ
    return environ.get(name, default)


# Paths #######################################################################

# The path temlate from which visual area labels are loaded. The template must
# match an MGZ or similar file that neuropythy can natively load that contains
# a vector of integer labels (values 0-10), one value per vertex.
# The template can contain the format strings {hemisphere} and {sid}.
label_path = getenv(
    "CMAG_LABEL_PATH",
    "/data/hcp/labels/visual/{hemisphere}.{sid}.mgz")

# The path in which to store cache data.
cache_path = getenv(
    "CMAG_CACHE_PATH",
    "/data/crowding/hcp-cmag")


# Labels ######################################################################

# The key for the labels in the label_path. The files loaded from the
# label_path must use these labels.
label_names = (
    # Label 0 is always the background label.
    None, 
    # Early visual cortex (1, 2, 3).
    'V1', 'V2', 'V3', 
    # Ventral visual cortex (4, 5, 6).
    'hV4', 'VO1', 'VO2',
    # Dorsal visual cortex (7, 8, 9, 10).
    'V3a', 'V3b', 'IPS0', 'LO1')
def label_lookup(name):
    """Returns the label number that matches the given label name.
    
    `label_lookup(name)` looks up the given label name (a string, for example
    `'V2'`) in the `cmag.config.label_names` sequence and returns its index. If
    the `label_names` entry in `cmag.config` is changed, this will continue to
    work.

    If the name is not found, then a `KeyError` is raised. Note that `None` is
    considered a valid name and will result in the return value of 0.
    """
    if label_names is not label_lookup._label_names:
        # Generate a new key.
        label_lookup._label_key = {
            (nm.upper() if isinstance(nm, str) else nm): index
            for (index, nm) in enumerate(label_names)}
        label_lookup._label_names = label_names
    if isinstance(name, str):
        name = name.upper()
    return label_lookup._label_key[name.upper()]


# Visual Field ################################################################

# The maximum eccentricity of the human field of view in degrees.
# This parameter is used to normalize various calculations related to models of
# cortical magnification. Because the magnification is very low near the
# periphery in any model, it typically isn't very sensitive to small changes.
fov = 200

# The maximum eccentricity assumed measured in a visual area.
# This parameter is used to normalize the eccentricity found in a visual area.
max_eccen = 7
