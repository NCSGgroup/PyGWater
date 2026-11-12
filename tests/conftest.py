import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# PyGWater's code is the folder src/, imported as `pygwater`. Normally it is installed (`pip install -e .`); if not,
# register src/ under that name so that the tests also run from a plain checkout (e.g. on GitHub CI).
try:
    import pygwater  # noqa: F401
except ImportError:
    spec = importlib.util.spec_from_file_location('pygwater', ROOT / 'src' / '__init__.py',
                                                  submodule_search_locations=[str(ROOT / 'src')])
    module = importlib.util.module_from_spec(spec)
    sys.modules['pygwater'] = module
    spec.loader.exec_module(module)

_NEEDS_LLN = ('test_consistency.py', 'test_basis.py')


def pytest_collection_modifyitems(config, items):
    """Tests that read the Love numbers are skipped when <data root>/input/LLN is missing (e.g. on GitHub CI)."""
    from pygwater.config import FileTool
    try:
        if (FileTool.get_input_dir('LLN') / 'PREM-LLNs-complete.dat').exists():
            return
    except (FileNotFoundError, ValueError):      # data root not set (e.g. GitHub CI)
        pass
    skip = pytest.mark.skip(reason='needs the load Love numbers in <data root>/input/LLN (see doc/installation.md)')
    for item in items:
        if item.fspath.basename in _NEEDS_LLN:
            item.add_marker(skip)
