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

import astro_metadata_translator


class BasicTestCase(unittest.TestCase):

    def test_basic(self):
        version = astro_metadata_translator.__version__
        self.assertIsNotNone(version)


if __name__ == "__main__":
    unittest.main()
