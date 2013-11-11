"""

Extension modules for fapswitch that do not add core functionality.

Add features that are not essential or require some external libraries that
are not standard.

"""

from fapswitch.config import warning

try:
    from .synthetic_accessibility import sa_score
except ImportError:
    warning("Unable to import SA Score; all scores will be zero.")
    def sa_score(*args):
        "Dummy SA Score that always gives 0"
        return 0.0
