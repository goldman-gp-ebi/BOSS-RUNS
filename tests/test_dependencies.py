import pytest
import logging
import stat
from pathlib import Path

from boss.dependencies import Dependencies, find_exe


def test_check_dependencies():
    dep = Dependencies()
    logging.info(dep.dependencies)
    logging.info(dep.minimap2)
    assert dep.minimap2 is not None
    # test if executable for user
    assert Path(dep.minimap2).stat().st_mode & stat.S_IXUSR



# Test find exe with different strings and expectations of finding an executable
find_exe_params = pytest.mark.parametrize("name, exp", [
    ("miniasm", True),
    ("dummy", False)
])

@find_exe_params
def test_find_exe(name, exp):
    exe = find_exe(name)
    assert bool(exe) is exp


