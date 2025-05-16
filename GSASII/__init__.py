import os
import sys

sys.path.append('.')

class PathHackingException(Exception): pass

def make_path_watcher():
    import sys
    import traceback

    init_path = tuple(sys.path)
    def no_path_hacking(event_name, args):
        if event_name != 'import':
            return
        module, filename, path, meta_path, path_hooks = args
        if path is not None and tuple(path) != init_path:
            lines = list(traceback.format_stack())
            for line in lines:
                print(line.strip())
            print(set(path) - set(init_path))
            #sys.exit(1)
            raise PathHackingException(lines[-2].strip())
    return no_path_hacking

if os.environ.get("GSASII_NOPATHHACKING", ''):
    sys.addaudithook(make_path_watcher())

del sys, make_path_watcher, os


from .GSASIIGUI import main  # noqa: E402

__all__ = ['main']
