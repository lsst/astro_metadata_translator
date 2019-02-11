"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import astro_metadata_translator
import astro_metadata_translator.version


globals().update(build_package_configs(
    project_name="astro_metadata_translator",
    version=astro_metadata_translator.version.__version__))


# Remove astro_metadata_translator from the default intersphinx configuration
try:
    del intersphinx_mapping['astro_metadata_translator']  # noqa
except KeyError:
    pass


# Add pipelines.lsst.io to the intersphinx configuration.
# NOTE: we might want to be more sophisticated about mapping corresponding
# versions of the Pipelines and astro_metadata_translator. This technique will
# mostly work if the Pipelines and astro_metadata_translator are developed
# concurrently.
intersphinx_mapping['lsst'] = ('https://pipelines.lsst.io/v/daily/', None)  # noqa
