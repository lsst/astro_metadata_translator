v30.0.4 (2026-02-28)
====================

New Features
------------

- Switched ``ObservationInfo`` implementation so that it is now a Pydantic ``BaseModel``.
  There should be no change from a user perspective. (`DM-54087 <https://rubinobs.atlassian.net/browse/DM-54087>`_)

- * Added new ``VisitInfoTranslator``, designed to be able to extract a minimal consistent set of information from the headers of visit images that were created from an ``lsst.afw.image.VisitInfo``.
    This translator is not enabled by default since it can result in confusion when an instrument translator can also translate the bulk of the headers.
  * Added new ``--translator-name`` command-line option to ``astrometadata translate`` to force a specific translator to be used (such as ``VisitInfo``).
    This does require that the translator has been registered (such as with the ``-p`` option).
  * Added ``MetadataTranslator.get_translator_by_name()`` class method to return a translator class from the name of the translator. (`DM-54255 <https://rubinobs.atlassian.net/browse/DM-54255>`_)

API Changes
-----------

- Added a new ``quiet`` parameter to ``ObservationInfo`` constructor and ``ObservationInfo.from_header()`` constructor.
  This can be used to turn warning log messages on translator failure into debug log messages. (`DM-54279 <https://rubinobs.atlassian.net/browse/DM-54279>`_)

Other Changes and Additions
---------------------------

- Migrated the documentation build to use ``sphinxutils``. (`DM-54087 <https://rubinobs.atlassian.net/browse/DM-54087>`_)


v30.0.0 (2026-01-15)
====================

Dropped support for Python 3.10.

New features
------------

* Modified the file readers (including command-line tooling) to support URIs as well as local files.
  This can be used to read headers from S3 buckets or HTTPS servers.  (`DM-41256 <https://rubinobs.atlassian.net/browse/DM-41256>`_)
* Added two new translated properties:

  - ``exposure_time_requested`` is the requested exposure time and ``exposure_time`` is now defined to be the actual exposure time.
    The default value for an older translator is for the two values to be the same.
  - ``altaz_end`` now contains the telescope position at the end of the observation.
    If a translator has not been updated to support this property the default value is `None`.
    (`DM-50382 <https://rubinobs.atlassian.net/browse/DM-50382>`_)

Miscellaneous Changes
---------------------

* Added a default clipping of the Alt/Az values, which can be overridden by translators.  (`DM-50167 <https://rubinobs.atlassian.net/browse/DM-50167>`_)


v29.0.0 (2025-04-16)
====================

Dropped support for Python 3.9.

New features
------------

* Added support for entry points to be declared by packages providing their own metadata translators.
  Use ``astrometadata --list-plugins`` to list all registered translators.
  Given the high overheads that can be encountered during the importing of translators, plugins have to be requested explicitly by label rather than preemptively loading all plugins.
  (`DM-47972 <https://rubinobs.atlassian.net/browse/DM-47972>`_)


Miscellaneous Changes
---------------------

* Improved the robustness of table mode output when there is a translation failure.
  (`DM-46970 <https://rubinobs.atlassian.net/browse/DM-46970>`_)


v28.0.0 (2025-01-24)
====================

Miscellaneous Changes
---------------------

* Improved robustness of table output when there are translation failures.
  (`DM-46797 <https://rubinobs.atlassian.net/browse/DM-46797>`_)
* Improved some DECam translations.
  (`DM-46014 <https://rubinobs.atlassian.net/browse/DM-46014>`_)
* Improved robustness of tests when the network connection fails.
  (`DM-44651 <https://rubinobs.atlassian.net/browse/DM-44651>`_)

v27.0.0 (2024-06-25)
====================

New features
------------

* Added new property, ``can_see_sky`` to report whether the observation can see sky photons or not.
  The default implementation tries to infer the answer from the observation types.
  (`DM-43103 <https://rubinobs.atlassian.net/browse/DM-43103>`_)

Miscellaneous Changes
---------------------

* Added a helper utility for calculating the observing day.
  (`DM-43109 <https://rubinobs.atlassian.net/browse/DM-43109>`_)
* Added support for ``observing_day`` and ``observing_day_offset`` properties.
  The default implementaiton assumes that the observing day is the UTC day.
  (`DM-42636 <https://rubinobs.atlassian.net/browse/DM-42636>`_)

v26.0.0 (2023-12-14)
====================

Drop support for Python 3.8.

Miscellaneous Changes
---------------------

* The baseline FITS metadata translator has been modified to prefer ``DATE-BEG`` over ``DATE-OBS``.
  (`DM-39731 <https://rubinobs.atlassian.net/browse/DM-39731>`_)
* Modified the ``assertObservationInfoFromYaml`` test helper method to support the ability to check Alt/Az values.
  (`DM-38400 <https://rubinobs.atlassian.net/browse/DM-38400>`_)
