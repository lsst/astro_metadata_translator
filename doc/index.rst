.. _astro_metadata_translator-package:

.. Title is the EUPS package name

#########################
astro_metadata_translator
#########################

.. Sentence/short paragraph describing what the package is for.

The ``astro_metadata_translator`` package provides generalized infrastructure
for handling metadata extraction for astronomical instrumentation.

There are header translation classes implemented as subclasses of
`~astro_metadata_translator.MetadataTranslator`.  These translation subclasses
implement methods corresponding to each derived property defined in
`~astro_metadata_translator.ObservationInfo`. The methods are named
`to_{property}` and can be implemented explicitly by a translation class, or
implicitly by defining trivial mappings from a header item to a property, or
constant mappings that are fixed for all headers independent of any header
values.
Defining a new translator subclass that inherits from
`~astro_metadata_translator.MetadataTranslator` and giving it a name,
automatically registers the translator as being available for automated header
translation.
A translation class does not need to reside in the
``astro_metadata_translator`` package.

`~astro_metadata_translator.ObservationInfo` is a class summarizing the
information from the translators.
An instance of this class can be instantiated from
any `dict`-like header.  By default the header translation class to use
is determined by asking each registered translator whether it knows how to
translate it.  If an explicit translation class should be used it can be
specified explicitly.

For details on the format of header correction files see `~astro_metadata_translator.fix_header`.

.. warning::
  The existing set of property names in
  `~astro_metadata_translator.ObservationInfo` should be considered as beta
  quality.
  Some of the names could yet be changed for consistency with other data
  dictionaries.


Project info
============

Repository
   https://github.com/lsst/astro_metadata_translator

JIRA component
   https://jira.lsstcorp.org/issues/?jql=component%3Dastro_metadata_translator

Command Line Utilities
======================

.. click:: astro_metadata_translator.cli.astrometadata:main
  :prog: astrometadata
  :show-nested:

Python API reference
====================

.. automodapi:: astro_metadata_translator

.. automodapi:: astro_metadata_translator.bin.translateheader

.. automodapi:: astro_metadata_translator.indexing
