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

"""Metadata translation code for afw VisitInfo headers."""

from __future__ import annotations

__all__ = ("VisitInfoTranslator",)

import logging
from collections.abc import Mapping, MutableMapping
from typing import TYPE_CHECKING, Any

import astropy.time
import astropy.units as u
from astropy.coordinates import EarthLocation

from ..translator import cache_translation
from .fits import FitsTranslator
from .helpers import altaz_from_degree_headers, tracking_from_degree_headers

if TYPE_CHECKING:
    import astropy.coordinates

log = logging.getLogger(__name__)


class VisitInfoTranslator(FitsTranslator):
    """Metadata translator for afw VisitInfo serialization.

    VisitInfo is only defined for on-sky observations.

    VisitInfo does not encode the following properties:

    * observing_day
    * detector_serial
    * detector_unique_name
    * detector_group
    * detector_name

    Some values can be found if there is butler provenance:

    * detector_num
    * visit_id
    """

    name = "VisitInfo"
    """Name of this translation class"""

    supported_instrument = None
    """Does not correspond to a single instrument."""

    default_resource_root = None
    """Corrections are not supported by this translator."""

    _const_map = {
        "detector_group": None,
        "detector_unique_name": None,
        "detector_serial": None,
        "detector_exposure_id": None,
        "detector_name": None,
        "physical_filter": None,
    }

    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = {
        "exposure_time": ("EXPTIME", {"unit": u.s}),
        "dark_time": ("DARKTIME", {"unit": u.s}),
        "boresight_airmass": "BORE-AIRMASS",
        "boresight_rotation_angle": ("BORE-ROTANG", {"unit": u.deg}),
        "observation_id": "IDNUM",
        "exposure_id": "IDNUM",
        "object": "OBJECT",
        "science_program": "PROGRAM",
        "telescope": ("TELESCOP", {"default": "Unknown"}),
        "instrument": "INSTRUMENT",
        "relative_humidity": "HUMIDITY",
        "temperature": ("AIRTEMP", {"unit": u.deg_C}),
        "pressure": ("AIRPRESS", {"unit": u.Pa}),
        "has_simulated_content": "HAS-SIMULATED-CONTENT",
        "observation_type": "OBSTYPE",
        "focus_z": ("FOCUSZ", {"unit": u.mm}),
        "observation_reason": "REASON",
        # Rely on butler provenance.
        "detector_num": "LSST BUTLER DATAID DETECTOR",
    }

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        """Indicate whether this translation class can translate the
        supplied header.

        Always returns `False`. This translator has to be selected explicitly
        from context.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        return False

    @cache_translation
    def to_observation_counter(self) -> int:
        """Return the exposure ID as proxy for counter.

        Some exposure IDs encode a year and counter in the integer, others
        simply have an incrementing counter

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return self.to_exposure_id()

    @cache_translation
    def to_boresight_rotation_coord(self) -> str:
        # Docstring will be inherited.
        if self.is_key_ok("ROTTYPE"):
            self._used_these_cards("ROTTYPE")
            return self._header["ROTTYPE"].lower()
        return "unknown"

    @cache_translation
    def to_visit_id(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        # This can be found in butler provenance in some cases.
        # It is not generally used though and visit_id is effectively
        # deprecated.
        prov_key = "LSST BUTLER DATAID VISIT"
        if self.is_key_ok(prov_key):
            self._used_these_cards(prov_key)
            return self._header[prov_key]
        return self.to_exposure_id()

    @cache_translation
    def to_datetime_begin(self) -> astropy.time.Time | None:
        if self.is_key_ok("DATE-AVG"):
            date_avg = self._from_fits_date("DATE-AVG", scale="tai")
            self._used_these_cards("DATE-AVG")
            return date_avg - (self.to_exposure_time() / 2.0)
        return None

    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        datetime_end = self.to_datetime_begin()
        if datetime_end is not None:
            datetime_end = self.to_datetime_begin() + self.to_exposure_time()
        return datetime_end

    @cache_translation
    def to_location(self) -> astropy.coordinates.EarthLocation:
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        # OBS-LONG is east positive for VisitInfo.
        lon = self._header["OBS-LONG"]
        value = EarthLocation.from_geodetic(lon, self._header["OBS-LAT"], self._header["OBS-ELEV"])
        self._used_these_cards("OBS-LONG", "OBS-LAT", "OBS-ELEV")
        return value

    @cache_translation
    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord | None:
        # Docstring will be inherited. Property defined in properties.py
        radecpairs = (("BORE-RA", "BORE-DEC"),)
        return tracking_from_degree_headers(self, "ICRS", radecpairs, unit=(u.deg, u.deg))

    @cache_translation
    def to_altaz_begin(self) -> astropy.coordinates.AltAz | None:
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(self, (("BORE-ALT", "BORE-AZ"),), self.to_datetime_begin())

    @classmethod
    def fix_header(
        cls, header: MutableMapping[str, Any], instrument: str, obsid: str, filename: str | None = None
    ) -> bool:
        """Fix DECam headers.

        Parameters
        ----------
        header : `dict`
            The header to update.  Updates are in place.
        instrument : `str`
            The name of the instrument.
        obsid : `str`
            Unique observation identifier associated with this header.
            Will always be provided.
        filename : `str`, optional
            Filename associated with this header. May not be set since headers
            can be fixed independently of any filename being known.

        Returns
        -------
        modified : `bool`
            Returns `True` if the header was updated.

        Notes
        -----
        No fixes are applied for VisitInfo at this time.
        """
        modified = False
        return modified
