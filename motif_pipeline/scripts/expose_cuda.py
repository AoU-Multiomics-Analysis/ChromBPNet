"""Expose TensorFlow CUDA wheel libraries when Cromwell replaces ENTRYPOINT."""
import importlib.util
from pathlib import Path
import sys

site_packages = Path(importlib.util.find_spec('tensorflow').origin).parent.parent
libraries = list((site_packages / 'nvidia').glob('*/lib/*.so*'))
if not libraries:
    raise ValueError('TensorFlow CUDA wheel libraries are missing from the motif image')
for library in libraries:
    target = site_packages / 'tensorflow' / library.name
    if not target.exists():
        target.symlink_to(library)
for executable in (site_packages / 'nvidia').glob('*/bin/ptxas'):
    target = Path(sys.prefix) / 'bin/ptxas'
    if not target.exists():
        target.symlink_to(executable)

print(f'[cuda] Exposed {len(libraries)} TensorFlow CUDA wheel libraries', flush=True)
