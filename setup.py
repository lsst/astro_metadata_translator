#!/usr/bin/env python
from setuptools import setup

version = "0.1.0"
with open("./python/astro_metadata_translator/version.py", "w") as f:
    f.write(f"__version__='{version}'")

setup(
    name='astro_metadata_translator',
    version=f"{version}",
    description='A translator for astronomical metadata.',
    long_description="",
    author='AURA/LSST',
    author_email='tjenness@lsst.org',
    url='https://github.com/lsst/astro_metadata_translator',
    packages=['astro_metadata_translator'],
    package_dir={'': 'python'},
    license='GPLv3',
    install_requires=[
        'astropy'
    ],
    tests_require=['pytest', 'pytest-flake8',],
    scripts=["bin.src/translate_header.py"]
)
