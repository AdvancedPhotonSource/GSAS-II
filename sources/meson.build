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

# Fortan compilation:

#histogram2d_source = custom_target('histogram2dmodule.c',
#  input : ['histogram2d.for'],  # .f so no F90 wrappers
#  output : ['histogram2dmodule.c', 'histogram2d-f2pywrappers.f'],
#  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'histogram2d', '--lower']
#)

r = run_command(py, '-m', 'numpy.f2py', 'pypowder.for', '-m', 'pypowder', '--lower', check:true)
py.extension_module('pypowder',
  ['pypowder.for', 'pypowdermodule.c', 'pypowder-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'polymask.for', '-m', 'polymask', '--lower', check:true)
py.extension_module('polymask',
  ['polymask.for', 'polymaskmodule.c', 'polymask-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'pydiffax.for', '-m', 'pydiffax', '--lower', check:true)
py.extension_module('pydiffax',
  ['pydiffax.for', 'pydiffaxmodule.c', 'pydiffax-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'histogram2d.for', '-m', 'histogram2d', '--lower', check:true)
py.extension_module('histogram2d',
  ['histogram2d.for', 'histogram2dmodule.c', 'histogram2d-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'pack_f.for', '-m', 'pack_f', '--lower', check:true)
py.extension_module('pack_f',
  ['pack_f.for', 'pack_fmodule.c', 'pack_f-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'pytexture.for', '-m', 'pytexture', '--lower', check:true)
py.extension_module('pytexture',
  ['pytexture.for', 'pytexturemodule.c', 'pytexture-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'spotmask.for', '-m', 'spotmask', '--lower', check:true)
py.extension_module('spotmask',
  ['spotmask.for', 'spotmaskmodule.c', 'spotmask-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'pyspg.for', '-m', 'pyspg', '--lower', check:true)
py.extension_module('pyspg',
  ['pyspg.for', 'pyspgmodule.c', 'pyspg-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'unpack_cbf.for', '-m', 'unpack_cbf', '--lower', check:true)
py.extension_module('unpack_cbf',
  ['unpack_cbf.for', 'unpack_cbfmodule.c', 'unpack_cbf-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)

r = run_command(py, '-m', 'numpy.f2py', 'histosigma2d.for', '-m', 'histosigma2d', '--lower', check:true)
py.extension_module('histosigma2d',
  ['histosigma2d.for', 'histosigma2dmodule.c', 'histosigma2d-f2pywrappers.f'],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'GSASII'
)


subdir('k_vec_cython')
