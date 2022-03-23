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

import unittest
from astro_metadata_translator import ObservationInfo

"""Test that ObservationInfo can be extended simply"""


class MyObsInfo(ObservationInfo):
    PROPERTIES = dict(
        myInfo=("My special information", "str", str, None, None),
        **ObservationInfo.PROPERTIES
    )


class ExtensionTestCase(unittest.TestCase):
    def test_basic(self):
        """Test construction of extended ObservationInfo"""
        string = "something different"
        obsinfo = MyObsInfo.makeObservationInfo(myInfo=string)

        # Behaves like the original
        self.assertIsInstance(obsinfo, ObservationInfo)
        self.assertIsInstance(obsinfo, MyObsInfo)
        self.assertEqual(obsinfo.myInfo, string)

        # Docstrings have been set
        self.assertTrue(MyObsInfo.myInfo.__doc__.startswith(MyObsInfo.PROPERTIES["myInfo"][0]))
        self.assertTrue(MyObsInfo.telescope.__doc__.startswith(ObservationInfo.PROPERTIES["telescope"][0]))

        with self.assertRaises(TypeError):
            # Type checking is applied, like in the original
            MyObsInfo.makeObservationInfo(myInfo=12345)

        with self.assertRaises(KeyError):
            # The original should still have no knowledge of "myInfo"
            ObservationInfo.makeObservationInfo(myInfo=string)


if __name__ == "__main__":
    unittest.main()
