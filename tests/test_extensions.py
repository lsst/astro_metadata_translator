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

import pickle
import unittest

from astro_metadata_translator import ObservationInfo, PropertyDefinition, StubTranslator, makeObservationInfo

"""Test that extensions to the core set of properties works"""

FOO = "bar"
NUMBER = 12345


class DummyTranslator(StubTranslator):
    name = "dummy"
    supported_instrument = "dummy"
    extensions = dict(
        number=PropertyDefinition("A number", "int", int),
        foo=PropertyDefinition("A string", "str", str),
    )
    _const_map = {
        "ext_foo": FOO,
    }

    @classmethod
    def can_translate(cls, header, filename=None):
        return "INSTRUME" in header and header["INSTRUME"] == "dummy"

    def to_ext_number(self):
        """Return the combination on my luggage"""
        return NUMBER


class ExtensionsTestCase(unittest.TestCase):
    def setUp(self):
        self.header = dict(INSTRUME="dummy")
        # The test translator is incomplete so it will issue warnings
        # about missing translations. Catch them.
        with self.assertWarns(UserWarning):
            self.obsinfo = ObservationInfo(self.header)

    def assert_observation_info(self, obsinfo):
        """Check that the `ObservationInfo` is as expected"""
        self.assertIsInstance(obsinfo, ObservationInfo)
        self.assertEqual(obsinfo.ext_foo, FOO)
        self.assertEqual(obsinfo.ext_number, NUMBER)

    def test_basic(self):
        """Test construction of extended ObservationInfo"""
        # Behaves like the original
        self.assert_observation_info(self.obsinfo)

        copy = makeObservationInfo(extensions=DummyTranslator.extensions, ext_foo=FOO, ext_number=NUMBER)
        self.assertEqual(copy, self.obsinfo)

        with self.assertRaises(AttributeError):
            # Variable is read-only
            self.obsinfo.ext_foo = "something completely different"

        with self.assertRaises(KeyError):
            # Can't specify extension value without declaring extensions
            makeObservationInfo(foo="foobar")

        with self.assertRaises(TypeError):
            # Type checking is applied, like in the original
            makeObservationInfo(extensions=DummyTranslator.extensions, ext_foo=98765)

    def test_pickle(self):
        """Test that pickling works on ObservationInfo with extensions"""
        obsinfo = pickle.loads(pickle.dumps(self.obsinfo))
        self.assert_observation_info(obsinfo)

    def test_simple(self):
        """Test that simple representation works"""
        simple = self.obsinfo.to_simple()
        self.assertIn("ext_foo", simple)
        self.assertIn("ext_number", simple)
        obsinfo = ObservationInfo.from_simple(simple)
        self.assert_observation_info(obsinfo)

    def test_json(self):
        """Test that JSON representation works"""
        json = self.obsinfo.to_json()
        self.assertIn("ext_foo", json)
        self.assertIn("ext_number", json)
        obsinfo = ObservationInfo.from_json(json)
        self.assert_observation_info(obsinfo)


if __name__ == "__main__":
    unittest.main()
