"""

Core configuration and options.

Expose universal logging methods and setup
global configuration availability.

"""

__all__ = ['options']

from logging import critical, warning, error, debug, info

from .config import Options

options = Options()