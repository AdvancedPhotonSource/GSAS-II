project(
    'gsas_ii',
    'c', 'cpp', 'cython', 'fortran',
    version: '0.1.0',
    license: 'BSD',
    meson_version: '>= 1.1.0',
    default_options: [
        'cpp_std=c++11',
    ],
)

# Seek the backend
if meson.backend() != 'ninja'
  error('Ninja backend required')
endif

# define some compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cy = meson.get_compiler('cython')
fc = meson.get_compiler('fortran')

m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif


if fc.get_id() == 'gcc'
  add_global_arguments('-fbounds-check', language : 'fortran')
endif

args = fc.first_supported_argument(['-ffree-form', '-free', '/free'])
assert(args != [], 'No arguments found?')


py = import('python').find_installation(pure: false)
py_dep = py.dependency()
pybind11_dep = dependency('pybind11')
numpy_dep = dependency('numpy')

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)


subdir('GSASII')
subdir('sources')
