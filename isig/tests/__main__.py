#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
__main__

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2016
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lgpl.html)
"""

if __name__ == "__main__":
    import inspect
    import os
    import pytest
    import sys
    PATH = os.path.dirname(os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe()))))
    sys.exit(pytest.main(args=[PATH, "--mpl"]))
