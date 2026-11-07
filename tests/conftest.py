import sys
from pathlib import Path

# make the src_* packages importable when pytest is started from anywhere
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


import pytest  # noqa: E402

_NEEDS_LLN = ('test_consistency.py', 'test_basis.py')


def pytest_collection_modifyitems(config, items):
    """Tests that read the Love numbers are skipped when <data root>/LLN is missing (e.g. on GitHub CI)."""
    from src_auxiliary.FileTool import FileTool
    if (FileTool.get_data_dir('LLN') / 'PREM-LLNs-complete.dat').exists():
        return
    skip = pytest.mark.skip(reason='needs the load Love numbers in <data root>/LLN (see doc/installation.md)')
    for item in items:
        if item.fspath.basename in _NEEDS_LLN:
            item.add_marker(skip)
