"""Compatibility module for the historic :mod:`Finder_charts` import path.

Use ``finder_charts`` for all new code. This module will be removed in a future
major release.
"""

from warnings import warn

warn(
    "Importing 'Finder_charts' is deprecated; import 'finder_charts' instead.",
    DeprecationWarning,
    stacklevel=2,
)

from finder_charts import Crosshair, Survey, Target, create_finder_charts

__all__ = ["Crosshair", "Survey", "Target", "create_finder_charts"]
