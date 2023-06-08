"""
Run this module to do pylint check.
"""

import os
from pathlib import Path


def pylint_check() -> int:
    rvt_library_dir_path = Path(__file__).parent.parent / "rvt"
    exit_code = os.system(f"pylint {rvt_library_dir_path}")
    if exit_code != 0 or exit_code != 4:  # 4 means warning
        raise Exception


if __name__ == '__main__':
    pylint_check()
