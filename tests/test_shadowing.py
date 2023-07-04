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

from astro_metadata_translator import StubTranslator


class ShadowBase(StubTranslator):
    """Base class for testing shadowing."""

    def to_instrument(self):
        return "BaseInstrument"


class ConstTranslator(StubTranslator):
    """Simple translation class with a constant mapping."""

    _const_map = {"instrument": "InstrumentB"}


class TrivialTranslator(ConstTranslator):
    """Translator inheriting from the constant mapping class but with an
    override.

    This class should not pick up the _const_map from parent class.
    """

    _trivial_map = {"instrument": "INSTRUME"}


class ExplicitTranslator(TrivialTranslator):
    """Translation class with explicit method inheriting from class
    with automatic translation methods.

    The explicit method should override the parent implementations
    and not inherit the _trivial_map from parent.
    """

    def to_instrument(self):
        return "InstrumentE"


class TranslatorShadowing(unittest.TestCase):
    """Test that shadowed translations are detected."""

    def test_shadowing(self):
        with self.assertLogs("astro_metadata_translator", level="WARN") as cm:

            class ShadowTranslator(StubTranslator):
                _const_map = {"instrument": "InstrumentC"}
                _trivial_map = {"instrument": "INSTRUME"}

                def to_instrument(self):
                    return "Instrument3"

        self.assertIn("defined in both", cm.output[0])
        self.assertIn("replaced by _const_map", cm.output[1])

        s = ShadowTranslator({})
        self.assertEqual(s.to_instrument(), "InstrumentC")

        with self.assertLogs("astro_metadata_translator", level="WARN") as cm:

            class ShadowTranslator(StubTranslator):
                _trivial_map = {"instrument": "INSTRUME"}

                def to_instrument(self):
                    return "Instrument3"

        self.assertIn("replaced by _trivial_map", cm.output[0])

        s = ShadowTranslator({"INSTRUME": "InstrumentT"})
        self.assertEqual(s.to_instrument(), "InstrumentT")

    def test_auto_maps1(self):
        t = TrivialTranslator({"INSTRUME": "InstrumentX"})
        self.assertEqual(t.to_instrument(), "InstrumentX")

    def test_auto_maps2(self):
        t = ExplicitTranslator({})
        self.assertEqual(t.to_instrument(), "InstrumentE")


if __name__ == "__main__":
    unittest.main()
