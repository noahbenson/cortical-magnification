# -*- coding: utf-8 -*-
###############################################################################
# __init__.py

"""A library for calculating cortical magnification.
"""

from .hh91 import (HH91, HH91_integral, HH91_find_a)
from .fitting import (fit_cumarea)

from . import hcp
from . import models
