# This file is part of obs_metadata.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Extract standard metadata from instrument FITS headers"""

__all__ = ("ObservationInfo", )

import logging

from .translator import MetadataTranslator

log = logging.getLogger(__name__)


class ObservationInfo:
    """Standardized representation of an instrument FITS header for a single
    exposure observation.

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument FITS header accessible as a `dict`.
    translator_class : `MetadataTranslator`-class, `optional`
        If not `None`, the class to use to translate the supplied headers
        into standard form. Otherwise each registered translator class will
        be asked in turn if it knows how to translate the supplied header.

    Raises
    ------
    ValueError
        The supplied header was not recognized by any of the registered
        translators.
    TypeError
        The supplied translator class was not a MetadataTranslator.
    """

    _PROPERTIES = {"telescope": ("Full name of the telescope", "str"),
                   "instrument": ("The instrument used to observe the Exposure", "str"),
                   "exposure": ("Unique (with instrument) integer identifier for this Exposure", "int"),
                   "visit": ("""ID of the Visit this Exposure is associated with.

    Science observations should essentially always be
    associated with a visit, but calibration observations
    may not be.""", "int"),
                   "abstract_filter": ("Generic name of this filter", "str"),
                   "physical_filter": ("The bandpass filter used for all exposures in this Visit.", "str"),
                   "snap": ("If visit is not null, the index of this Exposure in the Visit,"
                            "starting from zero.", "int"),
                   "datetime_begin": ("Timestamp of the start of the Exposure.", "astropy.time.Time"),
                   "datetime_end": ("Timestamp of the end of the Exposure.", "astropy.time.Time"),
                   "exposure_time": ("Duration of the Exposure with shutter open (seconds).", "float"),
                   "dark_time": ("Duration of the Exposure with shutter closed (seconds).", "float"),
                   "obsid": ("Unique observation identifier", "str")}

    """All the properties supported by this class with associated
    documentation."""

    def __init__(self, header, translator_class=None):
        if translator_class is None:
            translator_class = MetadataTranslator.determineTranslator(header)
        elif not issubclass(translator_class, MetadataTranslator):
            raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

        # Create an instance for this header
        translator = translator_class(header)

        # Loop over each translation (not final form -- this should be
        # defined in one place and consistent with translation classes)

        for t in self._PROPERTIES:
            # prototype code
            method = f"to_{t}"
            property = f"_{t}"

            try:
                print(f"Assigning {property} via {method}")
                setattr(self, property, getattr(translator, method)())
            except (AttributeError, KeyError):
                # For now assign None
                setattr(self, property, None)


# Dynamically add the standard properties
def _makeProperty(property, doc, return_type):
    def getter(self):
        return getattr(self, f"_{property}")

    getter.__doc__ = f"""{doc}

    Returns
    -------
    {property} : `{return_type}`
        Access the property.
    """
    return getter


for p, d in ObservationInfo._PROPERTIES.items():
    setattr(ObservationInfo, f"_{p}", None)
    setattr(ObservationInfo, p, property(_makeProperty(p, *d)))
