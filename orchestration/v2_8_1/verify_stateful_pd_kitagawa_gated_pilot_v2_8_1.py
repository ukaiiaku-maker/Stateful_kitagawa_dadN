#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path
here = Path(__file__).resolve().parent
raise SystemExit(subprocess.call([sys.executable, "-B", str(here / "test_v2_8_1_summary_discovery.py")]))
