"""
detettore6110: Detect and characterize insertion sequence polymorphisms and copy numbers in bacterial genomes.
"""

__version__ = "0.1.0"
__author__ = "Cristobal Strittmatter"
__license__ = "GPL-3.0"

from detettore6110.find import main as find_main
from detettore6110.summarize import main as summarize_main
from detettore6110.entry_point import main as cli_main

__all__ = [
    "__version__",
    "__author__",
    "__license__",
    "find_main",
    "summarize_main",
    "cli_main",
]
