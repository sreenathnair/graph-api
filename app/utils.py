#!/usr/bin/env python

"""
utils.py: Contains some useful utility functions
"""

import more_itertools as mit

__author__          = "Sreenath Sasidharan Nair"
__version__         = "1.0"
__email__           = "sreenath@ebi.ac.uk"

def find_ranges(iterable):
    """Yield range of consecutive numbers."""
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1]

            