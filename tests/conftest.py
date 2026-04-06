#!/usr/bin/env python3
"""Generate mock data files needed for tests."""
import subprocess
import sys
from pathlib import Path

# Generate mock FASTQs
gen = Path(__file__).parent / "data" / "generate_mock_data.py"
subprocess.run([sys.executable, str(gen)], check=True)
