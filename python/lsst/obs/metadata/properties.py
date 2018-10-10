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

"""Properties calculated by this package.

Defines all properties in one place so that both `ObservationInfo` and
`MetadataTranslator` can use them.  In particular, the translator
base class can use knowledge of these properties to predefine translation
stubs with documentation attached, and `ObservationInfo` can automatically
define the getter methods.

"""

__all__ = ("PROPERTIES", )


PROPERTIES = {"telescope": ("Full name of the telescope.", "str"),
              "instrument": ("The instrument used to observe the exposure.", "str"),
              "location": ("Location of the observatory.", "astropy.coordinates.EarthLocation"),
              "exposure": ("Unique (with instrument) integer identifier for this observation.", "int"),
              "visit": ("""ID of the Visit this Exposure is associated with.

Science observations should essentially always be
associated with a visit, but calibration observations
may not be.""", "int"),
              "abstract_filter": ("Generic name of this filter.", "str"),
              "physical_filter": ("The bandpass filter used for all exposures in this Visit.", "str"),
              "datetime_begin": ("Timestamp of the start of the Exposure.", "astropy.time.Time"),
              "datetime_end": ("Timestamp of the end of the Exposure.", "astropy.time.Time"),
              "exposure_time": ("Duration of the Exposure with shutter open (seconds).", "float"),
              "dark_time": ("Duration of the Exposure with shutter closed (seconds).", "float"),
              "boresight_airmass": ("Airmass of the boresight of the telescope.", "float"),
              "boresight_rotation_angle": ("Angle of the instrument in boresight_rotation_coord frame.",
                                           "astropy.coordinates.Angle"),
              "boresight_rotation_coord": ("Coordinate frame of the instrument rotation angle"
                                           " (options: sky, unknown).", "str"),
              "detector_num": ("Unique (for instrument) integer identifier for the sensor.", "int"),
              "detector_name": ("Name of the detector within the instrument (might not be unique).",
                                "str"),
              "detector_exposure_id": ("Unique integer identifier for this detector in this exposure.",
                                       "int"),
              "object": ("Object of interest or field name.", "str"),
              "temperature": ("Temperature outside the dome.", "astropy.units.Quantity"),
              "pressure": ("Atmospheric pressure outside the dome.", "astropy.units.Quantity"),
              "relative_humidity": ("Relative humidity outside the dome.", "float"),
              "tracking_radec": ("Requested RA/Dec to track.", "astropy.coordinates.SkyCoord"),
              "altaz_begin": ("Telescope boresight azimuth and elevation at start of observation.",
                              "astropy.coordinates.AltAz"),
              "science_program": ("Observing program (survey or proposal) identifier.", "str"),
              "obstype": ("Type of observation (currently: science, dark, flat, bias, focus).",
                          "str"),
              "obsid": ("Label uniquely identifying this observation (can be related to 'exposure').",
                        "str")}
