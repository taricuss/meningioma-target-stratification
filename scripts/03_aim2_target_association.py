#!/usr/bin/env python
# Aim 2: Target-gene program association
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from meningeal_extension.cli import aim2_main
if __name__ == "__main__":
    sys.exit(aim2_main())
