# This file is part of astro_metadata_translator.
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

"""Metadata translation code for CFHT MegaPrime FITS headers"""

__all__ = ("MegaPrimeTranslator", )

from astropy.coordinates import EarthLocation, SkyCoord, AltAz, Angle
import astropy.units as u

from .fits import FitsTranslator

filters = {'u.MP9301': 'u',
           'u.MP9302': 'u2',
           'g.MP9401': 'g',
           'g.MP9402': 'g2',
           'r.MP9601': 'r',
           'r.MP9602': 'r2',
           'i.MP9701': 'i',
           'i.MP9702': 'i2',
           'i.MP9703': 'i3',
           'z.MP9801': 'z',
           'z.MP9901': 'z2',
           }


class MegaPrimeTranslator(FitsTranslator):
    """Metadata translator for CFHT MegaPrime standard headers.
    """

    name = "MegaPrime"
    """Name of this translation class"""

    supported_instrument = "MegaPrime"
    """Supports the MegaPrime instrument."""

    _const_map = {"boresight_rotation_angle": Angle(float("nan")*u.deg),
                  "boresight_rotation_coord": "unknown"}

    _trivial_map = {"physical_filter": "FILTER",
                    "dark_time": ("DARKTIME", dict(unit=u.s)),
                    "exposure_time": ("EXPTIME", dict(unit=u.s)),
                    "observation_id": "OBSID",
                    "object": "OBJECT",
                    "science_program": "RUNID",
                    "exposure_id": "EXPNUM",
                    "visit_id": "EXPNUM",
                    "detector_name": "CCDNAME",
                    "relative_humidity": "RELHUMID",
                    "temperature": ("TEMPERAT", dict(unit=u.deg_C)),
                    "pressure": ("PRESSURE", dict(unit=u.hPa)),
                    "boresight_airmass": "AIRMASS"}

    def to_datetime_begin(self):
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            time_str=self._header["UTC-OBS"], scale="utc")
        self._used_these_cards("DATE-OBS", "UTC-OBS")
        return value

    def to_datetime_end(self):
        # Older files are missing UTCEND
        if "UTCEND" in self._header:
            # We know it is UTC
            value = self._from_fits_date_string(self._header["DATE-OBS"],
                                                time_str=self._header["UTCEND"], scale="utc")
            self._used_these_cards("DATE-OBS", "UTCEND")
        else:
            # Take a guess by adding on the exposure time
            value = self.to_datetime_begin() + self.to_exposure_time()
        return value

    def to_location(self):
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        # Height is not in MegaPrime files. Use the value from EarthLocation.of_site("CFHT")
        value = EarthLocation.from_geodetic(self._header["LONGITUD"], self._header["LATITUDE"], 4215.0)
        self._used_these_cards("LONGITUD", "LATITUDE")
        return value

    def to_detector_num(self):
        try:
            extname = self._header["EXTNAME"]
            num = int(extname[3:])  # chop off "ccd"
            self._used_these_cards("EXTNAME")
            return num
        except KeyError:
            # Dummy value, intended for PHU (need something to get filename)
            return 99

    def to_observation_type(self):
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        obstype = self._header["OBSTYPE"].strip().lower()
        self._used_these_cards("OBSTYPE")
        if obstype == "object":
            return "science"
        return obstype

    def to_tracking_radec(self):
        used = []
        for k in ("RADECSYS", "OBJRADEC", "RADESYS"):
            if k in self._header:
                frame = self._header[k].strip().lower()
                used.append(k)
                if frame == "gappt":
                    self._used_these_cards(*used)
                    # Moving target
                    return None
        else:
            frame = "icrs"
        radec = SkyCoord(self._header["RA_DEG"], self._header["DEC_DEG"],
                         frame=frame, unit=u.deg, obstime=self.to_datetime_begin(),
                         location=self.to_location())
        self._used_these_cards("RA_DEG", "DEC_DEG", *used)
        return radec

    def to_altaz_begin(self):
        altaz = AltAz(self._header["TELAZ"] * u.deg, self._header["TELALT"] * u.deg,
                      obstime=self.to_datetime_begin(), location=self.to_location())
        self._used_these_cards("TELALT", "TELAZ")
        return altaz

    def to_detector_exposure_id(self):
        return self.to_exposure_id() * 36 + self.to_detector_num()
