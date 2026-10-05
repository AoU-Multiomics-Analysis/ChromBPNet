"""Expose the pinned TensorFlow CUDA wheel libraries to the dynamic loader."""
import importlib.util
from pathlib import Path
import sys

site_packages = Path(importlib.util.find_spec('tensorflow').origin).parent.parent
for library in (site_packages / 'nvidia').glob('*/lib/*.so*'):
    target = site_packages / 'tensorflow' / library.name
    if not target.exists():
        target.symlink_to(library)
for executable in (site_packages / 'nvidia').glob('*/bin/ptxas'):
    target = Path(sys.prefix) / 'bin/ptxas'
    if not target.exists():
        target.symlink_to(executable)
