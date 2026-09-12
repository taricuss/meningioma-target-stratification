#!/usr/bin/env python
# Aim 4: Exploratory survival
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from meningeal_extension.cli import aim4_main
if __name__ == "__main__":
    sys.exit(aim4_main())
