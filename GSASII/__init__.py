import os
import sys

sys.path.append('.')


def make_path_watcher():
    import sys
    import traceback

    init_path = tuple(sys.path)
    def no_path_hacking(event_name, args):
        if event_name != 'import':
            return
        module, filename, path, meta_path, path_hooks = args
        if path is not None and tuple(path) != init_path:
            for line in traceback.format_stack():
                print(line.strip())
            print(set(path) - set(init_path))
            sys.exit(1)
    return no_path_hacking

if not os.environ.get("GSASII_YOLO_PATH", ''):
    sys.addaudithook(make_path_watcher())

del sys, make_path_watcher, os


from .GSASII import main  # noqa: E402

__all__ = ['main']
