py.extension_module('fmask',
        'fmask.c',
        dependencies : [py_dep, numpy_dep],
        install: true,
        subdir: 'GSASII',
        )

# Fortran compilation:
DIFFaXsubs_src = [
  'DIFFaXsubs/DIFFaXsubs.for',
]
DIFFaXsubs_lib = static_library('DIFFaXsubs_lib',
 DIFFaXsubs_src,
)
NISTlatsubs_src = [
  'NISTlatsubs/CENTER.f',
  'NISTlatsubs/CKPT02.f',
  'NISTlatsubs/DERIV.f',
  'NISTlatsubs/DETERM.f',
  'NISTlatsubs/DOT.f',
  'NISTlatsubs/HEAD0.f',
  'NISTlatsubs/HEAD2.f',
  'NISTlatsubs/HEAD3.f',
  'NISTlatsubs/HEAD4.f',
  'NISTlatsubs/HEAD5.f',
  'NISTlatsubs/HEAD6.f',
  'NISTlatsubs/HMATRX.f',
  'NISTlatsubs/INVERS.f',
  'NISTlatsubs/MNCOND.f',
  'NISTlatsubs/MULTIP.f',
  'NISTlatsubs/NORMAL.f',
  'NISTlatsubs/OUTPT1.f',
  'NISTlatsubs/OUTPT2.f',
  'NISTlatsubs/QMATRI.f',
  'NISTlatsubs/RD021.f',
  'NISTlatsubs/RD051.f',
  'NISTlatsubs/RD061.f',
  'NISTlatsubs/REDUCE.f',
  'NISTlatsubs/SET.f',
  'NISTlatsubs/SHORTV.f',
  'NISTlatsubs/SPCON2.f',
  'NISTlatsubs/SPCOND.f',
  'NISTlatsubs/TRANS.f',
  'NISTlatsubs/VOLUME.f',
]
NISTlatsubs_lib = static_library('NISTlatsubs_lib',
 NISTlatsubs_src,
)
powsubs_src = [
  'powsubs/acosd.for',
  'powsubs/cosd.for',
  'powsubs/epsvoigt.for',
  'powsubs/expint.for',
  'powsubs/gauleg.for',
  'powsubs/gerfc.for',
  'powsubs/hfunc.for',
  'powsubs/lorentz.for',
  'powsubs/psvfcj.f90',
  'powsubs/psvfcjexpo.for',
  'powsubs/psvfcjo.for',
  'powsubs/psvoigt.for',
  'powsubs/sind.for',
  'powsubs/tand.for',
]
powsubs_lib = static_library('powsubs_lib',
 powsubs_src,
)
texturesubs_src = [
  'texturesubs/aplms.for',
  'texturesubs/cosd.for',
  'texturesubs/dgammln.for',
  'texturesubs/factln.for',
  'texturesubs/plmpsi.for',
  'texturesubs/qlmn.for',
  'texturesubs/qlmninit.for',
  'texturesubs/sind.for',
]
texturesubs_lib = static_library('texturesubs_lib',
 texturesubs_src,
)


r = run_command(py, '-m', 'numpy.f2py', 'pypowder.for', '-m', 'pypowder', '--lower', check:true)
py.extension_module('pypowder',
  ['pypowder.for', 'pypowdermodule.c', 'pypowder-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
  link_with: powsubs_lib
)

r = run_command(py, '-m', 'numpy.f2py', 'pydiffax.for', '-m', 'pydiffax', '--lower', check:true)
py.extension_module('pydiffax',
  ['pydiffax.for', 'pydiffaxmodule.c', 'pydiffax-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
  link_with: DIFFaXsubs_lib
)

r = run_command(py, '-m', 'numpy.f2py', 'histogram2d.for', '-m', 'histogram2d', '--lower', check:true)
py.extension_module('histogram2d',
  ['histogram2d.for', 'histogram2dmodule.c', 'histogram2d-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
)

r = run_command(py, '-m', 'numpy.f2py', 'pack_f.for', '-m', 'pack_f', '--lower', check:true)
py.extension_module('pack_f',
  ['pack_f.for', 'pack_fmodule.c', 'pack_f-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
)

r = run_command(py, '-m', 'numpy.f2py', 'pytexture.for', '-m', 'pytexture', '--lower', check:true)
py.extension_module('pytexture',
  ['pytexture.for', 'pytexturemodule.c', 'pytexture-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
  link_with: texturesubs_lib
)

# create a copy of pyspg.for with the current git version info in FORTVER
pyspg_versioned = custom_target('pyspg_versioned',
  input : 'pyspg.for',
  output : 'pyspg_versioned.for',
  command : [py, files('patch_pyspg.py'), '@INPUT@', '@OUTPUT@'],
  build_always_stale : true,
)
r = run_command(py, '-m', 'numpy.f2py', 'pyspg.for', '-m', 'pyspg', '--lower', check:true)
py.extension_module('pyspg',
  [pyspg_versioned, 'pyspgmodule.c', 'pyspg-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
)

r = run_command(py, '-m', 'numpy.f2py', 'unpack_cbf.for', '-m', 'unpack_cbf', '--lower', check:true)
py.extension_module('unpack_cbf',
  ['unpack_cbf.for', 'unpack_cbfmodule.c', 'unpack_cbf-f2pywrappers.f'],
  dependencies : [py_dep, numpy_dep, fortranobject_dep],
  install : true,
  subdir: 'GSASII',
)

executable('convcell', files('convcell.f'), install:true)
executable('LATTIC', files('LATTIC.f','NISTlatsubs/BLOCKDATA1.f'), link_with:NISTlatsubs_lib, install:true)
subdir('k_vec_cython')

# create text file with version info
r = run_command(py, 'tagbinaries.py', check:true)

if host_machine.system() == 'windows'
   install = files('install.bat')
   sysinstall = files('sysinstall.bat')
else
   install = files('install.sh')
   sysinstall = files('sysinstall.sh')
endif 

r = run_target('local-install', command: install)
r = run_target('system-install', command: sysinstall)
