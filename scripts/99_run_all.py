#!/usr/bin/env python
# End-to-end pipeline: Aims 0-4 + figures + tables
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from meningeal_extension.cli import run_all
if __name__ == "__main__":
    sys.exit(run_all())
