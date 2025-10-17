"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.conf.pipelinespkg import *  # noqa: F403, import *

project = "astro_metadata_translator"
html_theme_options["logotext"] = project  # noqa: F405, unknown name
html_title = project
html_short_title = project
doxylink = {}

# Remove astro_metadata_translator from the default intersphinx configuration
try:
    del intersphinx_mapping["astro_metadata_translator"]  # noqa
except KeyError:
    pass


# Add pipelines.lsst.io to the intersphinx configuration.
# NOTE: we might want to be more sophisticated about mapping corresponding
# versions of the Pipelines and astro_metadata_translator. This technique will
# mostly work if the Pipelines and astro_metadata_translator are developed
# concurrently.
intersphinx_mapping["lsst"] = ("https://pipelines.lsst.io/v/daily/", None)  # noqa

nitpick_ignore_regex = [
    ("py:.*", r"lsst\..*"),  # Ignore warnings from links to other lsst packages.
    ("py:class", "(None|item) -- remove.*"),  # MutableSequence
    ("py:class", "integer -- return.*"),  # Sequence
    ("py:class", "ResourcePathExpression"),
    ("py:class", "ResourcePath"),
]
