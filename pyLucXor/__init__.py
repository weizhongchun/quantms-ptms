"""
pyLucXor - Python implementation of lucXor using PyOpenMS for mass spectrometry functionality.

This package provides tools for analyzing post-translational modifications (PTMs) in mass spectrometry data.
"""

__version__ = "0.1.0"

def __getattr__(name):
    if name == "LucXor":
        from .core import LucXor
        return LucXor
    elif name == "PyLuciPHOr2":
        from .cli import PyLuciPHOr2
        return PyLuciPHOr2
    elif name == "Peptide":
        from .peptide import Peptide
        return Peptide
    elif name == "PSM":
        from .psm import PSM
        return PSM
    elif name in ["CIDModel", "HCDModel"]:
        from .models import CIDModel, HCDModel
        return locals()[name]
    raise AttributeError(f"module 'pyLucXor' has no attribute '{name}'")

__all__ = ["LucXor", "PyLuciPHOr2", "Peptide", "PSM", "CIDModel", "HCDModel"]