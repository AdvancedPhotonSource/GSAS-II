project(
    'gsas_ii',
    'cpp',
    version: '0.1.0',
    license: 'BSD',
    meson_version: '>= 1.1.0',
    default_options: [
        'cpp_std=c++11',
    ],
)

py = import('python').find_installation(pure: false)
pybind11_dep = dependency('pybind11')

py.extension_module('_core',
    'src/main.cpp',
    subdir: 'gsas_ii',
    install: true,
    dependencies : [pybind11_dep],
)

install_subdir('src/gsas_ii', install_dir: py.get_install_dir() / 'gsas_ii', strip_directory: true)
