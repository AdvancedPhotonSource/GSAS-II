foreach file: ['DefaultExpressions.txt', 'OriginTemplate2.otpu', 'versioninfo.txt', 'Xsect.dat']

  install_data(
    file,
    install_tag: 'data',
    install_dir: py.get_install_dir(subdir: 'GSASII/inputs')
  )

endforeach

install_subdir(
  'GSASIImacros',
  install_tag: 'data',
  install_dir: py.get_install_dir(subdir: 'GSASII/inputs')
)
