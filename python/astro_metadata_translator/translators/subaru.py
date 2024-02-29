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

"""Metadata translation code for Subaru telescope."""

from __future__ import annotations

__all__ = ("SubaruTranslator",)

import astropy.time
from astropy.coordinates import EarthLocation

from ..translator import cache_translation
from .fits import FitsTranslator


class SubaruTranslator(FitsTranslator):
    """Metadata translator for Subaru telescope headers."""

    _observing_day_offset = astropy.time.TimeDelta(0, format="sec", scale="tai")

    @cache_translation
    def to_location(self) -> EarthLocation:
        """Return the location of the Subaru telescope on Mauna Kea.

        Hardcodes the location and does not look at any headers.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        return EarthLocation.from_geodetic(-155.476667, 19.825556, 4139.0)

    @cache_translation
    def to_observation_counter(self) -> int:
        """Return the lifetime exposure number.

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return self.to_exposure_id()

    @classmethod
    def observing_date_to_offset(cls, observing_date: astropy.time.Time) -> astropy.time.TimeDelta | None:
        """Return the offset to use when calculating the observing day.

        Parameters
        ----------
        observing_date : `astropy.time.Time`
            The date of the observation. Unused.

        Returns
        -------
        offset : `astropy.time.TimeDelta`
            The offset to apply. The offset is always 0 seconds. In Hawaii
            UTC rollover is at 2pm local time.
        """
        return cls._observing_day_offset
