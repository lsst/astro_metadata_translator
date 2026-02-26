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

import os.path
import unittest

import astropy.units as u

from astro_metadata_translator.tests import MetadataAssertHelper
from astro_metadata_translator.translators import VisitInfoTranslator

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class VisitInfoTestCase(unittest.TestCase, MetadataAssertHelper):
    """Test VisitInfo translations."""

    datadir = os.path.join(TESTDIR, "data")

    def test_decam_translator(self) -> None:
        test_data = (
            (
                "visit_info.yaml",
                dict(
                    telescope="Simonyi Survey Telescope",
                    instrument="LSSTCam",
                    boresight_rotation_coord="sky",
                    dark_time=30.9427 * u.s,
                    detector_exposure_id=None,
                    detector_name=None,
                    detector_unique_name=None,
                    detector_group=None,
                    detector_num=85,
                    detector_serial=None,
                    exposure_id=2025052000177,
                    exposure_group="2025052000177",
                    exposure_time=30.0011086463928 * u.s,
                    exposure_time_requested=30.0011086463928 * u.s,
                    focus_z=-1.8695301477484 * u.mm,
                    group_counter_end=2025052000177,
                    group_counter_start=2025052000177,
                    has_simulated_content=False,
                    object="COSMOS",
                    observation_counter=2025052000177,
                    observation_id="2025052000177",
                    observation_type="science",
                    observation_reason="field_survey_science",
                    observing_day=20250521,
                    physical_filter=None,
                    pressure=74510.0 * u.Pa,
                    relative_humidity=21.3250007629395,
                    science_program="BLOCK-365",
                    temperature=12.1499996185303 * u.deg_C,
                    visit_id=2025052000177,
                    wcs_params=dict(max_sep=3.5),
                ),
            ),
        )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoFromYaml(
                    file, dir=self.datadir, translator_class=VisitInfoTranslator, **expected
                )


if __name__ == "__main__":
    unittest.main()
