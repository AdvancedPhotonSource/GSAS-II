# Where is this source file ?
#py.extension_module('_core',
#    'src/main.cpp',
#    subdir: 'gsas_ii',
#    install: true,
#    dependencies : [pybind11_dep],
#)

#install_subdir('src/gsas_ii', install_dir: py.get_install_dir() / 'gsas_ii', strip_directory: true)

py.extension_module('fmask',
        'fmask.c',
        dependencies : [py_dep, numpy_dep],
#        include_directories : ['include'],
        install: true,
        subdir: 'GSASII',
        )

