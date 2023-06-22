"""
Package for analysis of electrophysiology data in Python.

"""

try:
    import importlib.metadata
    __version__ = importlib.metadata.version("pyelectro")
except ImportError:
    import importlib_metadata
    __version__ = importlib_metadata.version("pyelectro")
