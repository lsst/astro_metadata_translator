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

"""Transform ObservationInfo into "standard" FITS headers."""

__all__ = ("info_to_fits", "dates_to_fits", "group_to_fits")


def dates_to_fits(date_begin, date_end):
    """Convert two dates into FITS form.

    Parameters
    ----------
    date_begin : `astropy.time.Time`
        Date representing the beginning of the observation.
    date_end  : `astropy.time.Time`
        Date representing the end of the observation.

    Returns
    -------
    cards : `dict` of `str` to `str` or `float`
        Header card keys and values following the FITS standard.
        If neither date is defined this may be empty.
    """
    cards = {}
    if date_begin is None and date_end is None:
        # no date headers can be written
        return cards

    cards["TIMESYS"] = "TAI"

    date_avg = None
    if date_begin is not None and date_end is not None:
        date_avg = date_begin + (date_end - date_begin)/2.0

    for fragment, date in (("OBS", date_begin), ("END", date_end), ("AVG", date_avg)):
        if date is not None:
            tai = date.tai
            cards[f"DATE-{fragment}"] = tai.isot
            cards[f"MJD-{fragment}"] = tai.mjd

    return cards


def info_to_fits(obs_info):
    """Convert an `ObservationInfo` to something suitable for writing
    to a FITS file.

    Parameters
    ----------
    obs_info : `ObservationInfo`
        Standardized observation information to transform to FITS headers.

    Returns
    -------
    cards : `dict` of `str` to (`int`, `float`, `str`, `bool`)
        FITS header keys and values in form understood by FITS.
    comments : `dict` of `str` to `str`
        Suitable comment string.  There will be at most one entry for each key
        in ``cards``.
    """

    cards = {}
    comments = {}

    if obs_info.instrument is not None:
        cards["INSTRUME"] = obs_info.instrument
        comments["INSTRUME"] = "Name of instrument"

    cards.update(dates_to_fits(obs_info.datetime_begin, obs_info.datetime_end))

    return cards, comments


def group_to_fits(obs_group):
    """Convert an `ObservationGroup` to something suitable for writing
    to a FITS file.

    Parameters
    ----------
    obs_group : `ObservationGroup`
        Collection of observation information to transform to a single
        FITS header.

    Returns
    -------
    cards : `dict` of `str` to (`int`, `float`, `str`, `bool`)
        FITS header keys and values in form understood by FITS.
    comments : `dict` of `str` to `str`
        Suitable comment string.  There will be at most one entry for each key
        in ``cards``.
    """

    cards = {}
    comments = {}

    oldest, newest = obs_group.extremes()

    instruments = obs_group.property_values("instrument")
    if len(instruments) == 1:
        cards["INSTRUME"] = list(instruments)[0]
        comments["INSTRUME"] = "Name of instrument"

    cards.update(dates_to_fits(oldest.datetime_begin, newest.datetime_end))

    return cards, comments
