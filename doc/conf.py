"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import astro_metadata_translator
import astro_metadata_translator.version


_g = globals()
_g.update(build_package_configs(
    project_name="astro_metadata_translator",
    version=astro_metadata_translator.version.__version__))
