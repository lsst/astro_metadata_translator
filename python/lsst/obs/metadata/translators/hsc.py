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

"""Metadata translation code for HSC FITS headers"""

__all__ = ("HscTranslator", )

import re

from .subaru import SubaruTranslator


class HscTranslator(SubaruTranslator):
    """Metadata translator for HSC standard headers.
    """

    name = "HSC"
    """Name of this translation class"""

    supportedInstrument = "HSC"
    """Supports the HSC instrument."""

    _constMap = {"instrument": "HSC"}
    """This translator only works for HSC."""

    _trivialMap = {"physical_filter": "FILTER01",
                   "obsid": "EXP-ID",
                   "object": "OBJECT",
                   "science_program": "PROP-ID",
                   "detector_num": "DET-ID",
                   "detector_name": "T_CCDSN",
                   "boresight_airmass": "AIRMASS",
                   "exposure_time": "EXPTIME",
                   "dark_time": "EXPTIME",  # Assume same as exposure time
                   }
    """One-to-one mappings"""

    @classmethod
    def canTranslate(cls, header):
        """Indicate whether this translation class can translate the
        supplied header.

        There is no ``INSTRUME`` header in an HSC file, so this method
        looks for HSC mentions in other headers.

        Parameters
        ----------
        header : `dict`-like
           Header to convert to standardized form.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        for k in ("EXP-ID", "FRAMEID"):
            if k in header:
                if header[k].startswith("HSC"):
                    return True
        return False

    def to_datetime_begin(self):
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            timeStr=self._header["UT"], scale="utc")
        self._used_these_cards("DATE-OBS", "UT")
        return value

    def to_datetime_end(self):
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            timeStr=self._header["UT-END"], scale="utc")
        self._used_these_cards("DATE-OBS", "UT-END")
        return value

    def to_abstract_filter(self):
        physical = self.to_physical_filter()
        if physical.startswith("HSC-"):
            return physical[4:]
        return None

    def to_visit(self):
        """Calculate unique visit integer for this observation

        Returns
        -------
        visit : `int`
            Integer uniquely identifying this visit.
        """
        expId = self._header["EXP-ID"].strip()
        m = re.search("^HSCE(\d{8})$", expId)  # 2016-06-14 and new scheme
        if m:
            self._used_these_cards("EXP-ID")
            return int(m.group(1))

        # Fallback to old scheme
        m = re.search("^HSC([A-Z])(\d{6})00$", expId)
        if not m:
            raise RuntimeError(f"Unable to interpret EXP-ID: {expId}")
        letter, visit = m.groups()
        visit = int(visit)
        if visit == 0:
            # Don't believe it
            frameId = self._header["FRAMEID"].strip()
            m = re.search("^HSC([A-Z])(\d{6})\d{2}$", frameId)
            if not m:
                raise RuntimeError(f"Unable to interpret FRAMEID: {frameId}")
            letter, visit = m.groups()
            visit = int(visit)
            if visit % 2:  # Odd?
                visit -= 1
        self._used_these_cards("EXP-ID", "FRAMEID")
        return visit + 1000000*(ord(letter) - ord("A"))

    def to_exposure(self):
        """Calculate the unique integer ID for this exposure.

        Assumed to be identical to the visit ID in this implementation.

        Returns
        -------
        exp : `int`
            Unique exposure identifier.
        """
        return self.to_visit()
