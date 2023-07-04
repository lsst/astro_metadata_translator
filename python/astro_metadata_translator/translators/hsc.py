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

"""Metadata translation code for HSC FITS headers."""

from __future__ import annotations

__all__ = ("HscTranslator",)

import logging
import posixpath
import re
from collections.abc import Mapping
from typing import Any

import astropy.units as u
from astropy.coordinates import Angle

from ..translator import CORRECTIONS_RESOURCE_ROOT, cache_translation
from .suprimecam import SuprimeCamTranslator

log = logging.getLogger(__name__)


class HscTranslator(SuprimeCamTranslator):
    """Metadata translator for HSC standard headers."""

    name = "HSC"
    """Name of this translation class"""

    supported_instrument = "HSC"
    """Supports the HSC instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "HSC")
    """Default resource path root to use to locate header correction files."""

    _const_map = {"instrument": "HSC", "boresight_rotation_coord": "sky"}
    """Hard wire HSC even though modern headers call it Hyper Suprime-Cam"""

    _trivial_map = {
        "detector_serial": "T_CCDSN",
    }
    """One-to-one mappings"""

    # Zero point for HSC dates: 2012-01-01  51544 -> 2000-01-01
    _DAY0 = 55927

    # CCD index mapping for commissioning run 2
    _CCD_MAP_COMMISSIONING_2 = {
        112: 106,
        107: 105,
        113: 107,
        115: 109,
        108: 110,
        114: 108,
    }

    _DETECTOR_NUM_TO_UNIQUE_NAME = [
        "1_53",
        "1_54",
        "1_55",
        "1_56",
        "1_42",
        "1_43",
        "1_44",
        "1_45",
        "1_46",
        "1_47",
        "1_36",
        "1_37",
        "1_38",
        "1_39",
        "1_40",
        "1_41",
        "0_30",
        "0_29",
        "0_28",
        "1_32",
        "1_33",
        "1_34",
        "0_27",
        "0_26",
        "0_25",
        "0_24",
        "1_00",
        "1_01",
        "1_02",
        "1_03",
        "0_23",
        "0_22",
        "0_21",
        "0_20",
        "1_04",
        "1_05",
        "1_06",
        "1_07",
        "0_19",
        "0_18",
        "0_17",
        "0_16",
        "1_08",
        "1_09",
        "1_10",
        "1_11",
        "0_15",
        "0_14",
        "0_13",
        "0_12",
        "1_12",
        "1_13",
        "1_14",
        "1_15",
        "0_11",
        "0_10",
        "0_09",
        "0_08",
        "1_16",
        "1_17",
        "1_18",
        "1_19",
        "0_07",
        "0_06",
        "0_05",
        "0_04",
        "1_20",
        "1_21",
        "1_22",
        "1_23",
        "0_03",
        "0_02",
        "0_01",
        "0_00",
        "1_24",
        "1_25",
        "1_26",
        "1_27",
        "0_34",
        "0_33",
        "0_32",
        "1_28",
        "1_29",
        "1_30",
        "0_41",
        "0_40",
        "0_39",
        "0_38",
        "0_37",
        "0_36",
        "0_47",
        "0_46",
        "0_45",
        "0_44",
        "0_43",
        "0_42",
        "0_56",
        "0_55",
        "0_54",
        "0_53",
        "0_31",
        "1_35",
        "0_35",
        "1_31",
        "1_48",
        "1_51",
        "1_52",
        "1_57",
        "0_57",
        "0_52",
        "0_51",
        "0_48",
    ]

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        """Indicate whether this translation class can translate the
        supplied header.

        There is no ``INSTRUME`` header in early HSC files, so this method
        looks for HSC mentions in other headers.  In more recent files the
        instrument is called "Hyper Suprime-Cam".

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
        if "INSTRUME" in header:
            return header["INSTRUME"] == "Hyper Suprime-Cam"

        for k in ("EXP-ID", "FRAMEID"):
            if cls.is_keyword_defined(header, k):
                if header[k].startswith("HSC"):
                    return True
        return False

    @cache_translation
    def to_exposure_id(self) -> int:
        """Calculate unique exposure integer for this observation.

        Returns
        -------
        visit : `int`
            Integer uniquely identifying this exposure.
        """
        exp_id = self._header["EXP-ID"].strip()
        m = re.search(r"^HSCE(\d{8})$", exp_id)  # 2016-06-14 and new scheme
        if m:
            self._used_these_cards("EXP-ID")
            return int(m.group(1))

        # Fallback to old scheme
        m = re.search(r"^HSC([A-Z])(\d{6})00$", exp_id)
        if not m:
            raise RuntimeError(f"{self._log_prefix}: Unable to interpret EXP-ID: {exp_id}")
        letter, visit = m.groups()
        visit = int(visit)
        if visit == 0:
            # Don't believe it
            frame_id = self._header["FRAMEID"].strip()
            m = re.search(r"^HSC([A-Z])(\d{6})\d{2}$", frame_id)
            if not m:
                raise RuntimeError(f"{self._log_prefix}: Unable to interpret FRAMEID: {frame_id}")
            letter, visit = m.groups()
            visit = int(visit)
            if visit % 2:  # Odd?
                visit -= 1
        self._used_these_cards("EXP-ID", "FRAMEID")
        return visit + 1000000 * (ord(letter) - ord("A"))

    @cache_translation
    def to_boresight_rotation_angle(self) -> Angle:
        # Docstring will be inherited. Property defined in properties.py
        # Rotation angle formula determined empirically from visual inspection
        # of HSC images.  See DM-9111.
        angle = Angle(270.0 * u.deg) - Angle(self.quantity_from_card("INST-PA", u.deg))
        angle = angle.wrap_at("360d")
        return angle

    @cache_translation
    def to_detector_num(self) -> int:
        """Calculate the detector number.

        Focus CCDs were numbered incorrectly in the readout software during
        commissioning run 2.  This method maps to the correct ones.

        Returns
        -------
        num : `int`
            Detector number.
        """
        ccd = super().to_detector_num()
        try:
            tjd = self._get_adjusted_mjd()
        except Exception:
            return ccd

        if tjd > 390 and tjd < 405:
            ccd = self._CCD_MAP_COMMISSIONING_2.get(ccd, ccd)

        return ccd

    @cache_translation
    def to_detector_exposure_id(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id() * 200 + self.to_detector_num()

    @cache_translation
    def to_detector_group(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        unique = self.to_detector_unique_name()
        return unique.split("_")[0]

    @cache_translation
    def to_detector_unique_name(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        # Mapping from number to unique name is defined solely in camera
        # geom files.
        # There is no header for it.
        num = self.to_detector_num()
        return self._DETECTOR_NUM_TO_UNIQUE_NAME[num]

    @cache_translation
    def to_detector_name(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        # Name is defined from unique name
        unique = self.to_detector_unique_name()
        return unique.split("_")[1]

    @cache_translation
    def to_focus_z(self) -> u.Quantity:
        # Docstring will be inherited. Property defined in properties.py
        foc_val = self._header["FOC-VAL"]
        return foc_val * u.mm
