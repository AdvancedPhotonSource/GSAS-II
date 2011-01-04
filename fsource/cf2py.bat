path=%path:MinGW32=MinGW64%
f2py -c -m --opt='-O1' %1 %1.for %2 --compiler=mingw32
path=%path:MinGW64=MinGW32%
move %1.pyd ..\binwin64-2.6
