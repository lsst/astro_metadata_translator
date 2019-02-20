#!/usr/bin/env python
from setuptools import setup

version = "0.1.0"
with open("./python/astro_metadata_translator/version.py", "w") as f:
    print(f"""
__all__ = ("__version__", )
__version__='{version}'""", file=f)

setup(
    version=f"{version}",
)
