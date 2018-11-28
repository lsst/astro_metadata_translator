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

import unittest

from astro_metadata_translator import merge_headers


class HeadersTestCase(unittest.TestCase):

    def setUp(self):
        # Define reference headers
        self.h1 = dict(
            ORIGIN="LSST",
            KEY0=0,
            KEY1=1,
            KEY2=3,
            KEY3=3.1415,
            KEY4="a",
        )
        self.h2 = dict(
            ORIGIN="LSST",
            KEY0="0",
            KEY2=4,
            KEY5=42
        )
        self.h3 = dict(
            ORIGIN="AUXTEL",
            KEY3=3.1415,
            KEY2=50,
            KEY5=42,
        )
        self.h4 = dict(
            KEY6="New",
            KEY1="Exists",
        )

    def test_fail(self):
        with self.assertRaises(ValueError):
            merge_headers([self.h1], mode="wrong")

        with self.assertRaises(ValueError):
            merge_headers([])

    def test_merging_overwrite(self):
        merged = merge_headers([self.h1, self.h2], mode="overwrite")

        expected = {
            "ORIGIN": self.h1["ORIGIN"],
            "KEY0": self.h2["KEY0"],
            "KEY1": self.h1["KEY1"],
            "KEY2": self.h2["KEY2"],
            "KEY3": self.h1["KEY3"],
            "KEY4": self.h1["KEY4"],
            "KEY5": self.h2["KEY5"],
        }
        self.assertEqual(merged, expected)

        merged = merge_headers([self.h1, self.h2, self.h3, self.h4],
                               mode="overwrite")

        expected = {
            "ORIGIN": self.h3["ORIGIN"],
            "KEY0": self.h2["KEY0"],
            "KEY1": self.h4["KEY1"],
            "KEY2": self.h3["KEY2"],
            "KEY3": self.h3["KEY3"],
            "KEY4": self.h1["KEY4"],
            "KEY5": self.h3["KEY5"],
            "KEY6": self.h4["KEY6"],
        }

        self.assertEqual(merged, expected)

    def test_merging_first(self):
        merged = merge_headers([self.h1, self.h2, self.h3, self.h4],
                               mode="first")

        expected = {
            "ORIGIN": self.h1["ORIGIN"],
            "KEY0": self.h1["KEY0"],
            "KEY1": self.h1["KEY1"],
            "KEY2": self.h1["KEY2"],
            "KEY3": self.h1["KEY3"],
            "KEY4": self.h1["KEY4"],
            "KEY5": self.h2["KEY5"],
            "KEY6": self.h4["KEY6"],
        }

        self.assertEqual(merged, expected)

    def test_merging_drop(self):
        merged = merge_headers([self.h1, self.h2, self.h3, self.h4],
                               mode="drop")

        expected = {
            "KEY3": self.h1["KEY3"],
            "KEY4": self.h1["KEY4"],
            "KEY5": self.h2["KEY5"],
            "KEY6": self.h4["KEY6"],
        }

        self.assertEqual(merged, expected)

    def test_merging_append(self):
        # Try with two headers first
        merged = merge_headers([self.h1, self.h2], mode="append")

        expected = {
            "ORIGIN": self.h1["ORIGIN"],
            "KEY0": [self.h1["KEY0"], self.h2["KEY0"]],
            "KEY1": self.h1["KEY1"],
            "KEY2": [self.h1["KEY2"], self.h2["KEY2"]],
            "KEY3": self.h1["KEY3"],
            "KEY4": self.h1["KEY4"],
            "KEY5": self.h2["KEY5"],
        }

        self.assertEqual(merged, expected)

        merged = merge_headers([self.h1, self.h2, self.h3, self.h4],
                               mode="append")

        expected = {
            "ORIGIN": [self.h1["ORIGIN"], self.h2["ORIGIN"], self.h3["ORIGIN"], None],
            "KEY0": [self.h1["KEY0"], self.h2["KEY0"], None, None],
            "KEY1": [self.h1["KEY1"], None, None, self.h4["KEY1"]],
            "KEY2": [self.h1["KEY2"], self.h2["KEY2"], self.h3["KEY2"], None],
            "KEY3": self.h3["KEY3"],
            "KEY4": self.h1["KEY4"],
            "KEY5": self.h3["KEY5"],
            "KEY6": self.h4["KEY6"],
        }

        self.assertEqual(merged, expected)


if __name__ == "__main__":
    unittest.main()
