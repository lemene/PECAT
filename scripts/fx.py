#!/usr/bin/env python3

import os
import sys


#Setting executable paths
fx_root = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, fx_root)
os.environ["PATH"] = fx_root + os.pathsep + os.environ["PATH"]

from fxtools.main import main
sys.exit(main())
