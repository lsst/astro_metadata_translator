# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file at the top-level directory of this distribution
# for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

import logging as _logging
from importlib.metadata import entry_points as _entry_points

from .headers import *
from .observationGroup import *
from .observationInfo import *
from .properties import *
from .translator import *
from .translators import *
from .version import *

# Load registered translators. This runs after all the other package
# code has been imported.
_LOG = _logging.getLogger(__name__)
_plugins = _entry_points(group="astro_metadata_translators")
for _p in _plugins:
    # Only need to load the entry point since the packages will
    # automatically register themselves via __init_subclass__.
    # The only requirement is that the entry point code forces the
    # translator code to be imported.
    try:
        _p.load()
    except Exception as _e:
        _LOG.warning(f"Plugin {_p} failed to load. Skipping: {_e!r}")
