#!/usr/bin/env python3
# -*- coding: utf-8 -*-
try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError
    
try:
    __version__ = version("GoMartini")
except PackageNotFoundError:
    # If the package is not installed, don't add __version__
    pass