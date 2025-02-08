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

from __future__ import annotations

__all__ = ("dates_to_fits", "group_to_fits", "info_to_fits")

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import astropy.time

    from ..observationGroup import ObservationGroup
    from ..observationInfo import ObservationInfo


def dates_to_fits(date_begin: astropy.time.Time, date_end: astropy.time.Time) -> dict[str, Any]:
    """Convert two dates into FITS form.

    Parameters
    ----------
    date_begin : `astropy.time.Time`
        Date representing the beginning of the observation.
    date_end : `astropy.time.Time`
        Date representing the end of the observation.

    Returns
    -------
    cards : `dict` of `str` to `str` or `float`
        Header card keys and values following the FITS standard.
        If neither date is defined this may be empty.
    """
    cards: dict[str, Any] = {}
    if date_begin is None and date_end is None:
        # no date headers can be written
        return cards

    cards["TIMESYS"] = "TAI"

    date_avg = None
    if date_begin is not None and date_end is not None:
        date_avg = date_begin + (date_end - date_begin) / 2.0

    for fragment, date in (("OBS", date_begin), ("BEG", date_begin), ("END", date_end), ("AVG", date_avg)):
        if date is not None:
            tai = date.tai
            cards[f"DATE-{fragment}"] = tai.isot
            cards[f"MJD-{fragment}"] = tai.mjd

    return cards


def info_to_fits(obs_info: ObservationInfo) -> tuple[dict[str, Any], dict[str, str]]:
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


def group_to_fits(obs_group: ObservationGroup) -> tuple[dict[str, Any], dict[str, str]]:
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
