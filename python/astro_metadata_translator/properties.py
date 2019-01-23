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
              "exposure_id": ("Unique (with instrument) integer identifier for this observation.", "int"),
              "visit_id": ("""ID of the Visit this Exposure is associated with.

Science observations should essentially always be
associated with a visit, but calibration observations
may not be.""", "int"),
              "physical_filter": ("The bandpass filter used for this observation.", "str"),
              "datetime_begin": ("Time of the start of the observation.", "astropy.time.Time"),
              "datetime_end": ("Time of the end of the observation.", "astropy.time.Time"),
              "exposure_time": ("Duration of the exposure with shutter open (seconds).",
                                "astropy.units.Quantity"),
              "dark_time": ("Duration of the exposure with shutter closed (seconds).",
                            "astropy.units.Quantity"),
              "boresight_airmass": ("Airmass of the boresight of the telescope.", "float"),
              "boresight_rotation_angle": ("Angle of the instrument in boresight_rotation_coord frame.",
                                           "astropy.coordinates.Angle"),
              "boresight_rotation_coord": ("Coordinate frame of the instrument rotation angle"
                                           " (options: sky, unknown).", "str"),
              "detector_num": ("Unique (for instrument) integer identifier for the sensor.", "int"),
              "detector_name": ("Name of the detector within the instrument (might not be unique).",
                                "str"),
              "detector_serial": ("Serial number/string associated with this detector.", "str"),
              "detector_group": ("Collection name of which this detector is a part. "
                                 "Can be None if there are no detector groupings.", "str"),
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
              "observation_type": ("Type of observation (currently: science, dark, flat, bias, focus).",
                                   "str"),
              "observation_id": ("Label uniquely identifying this observation "
                                 "(can be related to 'exposure_id').",
                                 "str")}
