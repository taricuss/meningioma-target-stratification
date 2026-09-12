#!/usr/bin/env python
# Aim 0: Cohort lock-in
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from meningeal_extension.cli import aim0_main
if __name__ == "__main__":
    sys.exit(aim0_main())
