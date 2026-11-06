import sys
from pathlib import Path

# make the src_* packages importable when pytest is started from anywhere
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
