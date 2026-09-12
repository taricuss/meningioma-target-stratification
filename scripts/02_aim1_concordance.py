#!/usr/bin/env python
# Aim 1: Classifier concordance
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from meningeal_extension.cli import aim1_main
if __name__ == "__main__":
    sys.exit(aim1_main())
