import glob
import os
import importlib.util
import sys

__all__ = []
for f in sorted(glob.glob(os.path.join(os.path.dirname(__file__),'G2*.py'))):
    nam = os.path.splitext(os.path.split(f)[1])[0]
    try:
        modspec = importlib.util.spec_from_file_location(nam, f)
        module = importlib.util.module_from_spec(modspec)
        module.__package__ = "GSASII.imports"
        modspec.loader.exec_module(module)
        sys.modules[f"GSASII.imports.{nam}"] = module
        __all__.append(nam)
    except ImportError as msg:
        print(f'Failed to import importer {nam} with error\n{msg}')
