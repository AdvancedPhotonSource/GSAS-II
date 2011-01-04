f2py -c -m %1 %1.for %2 --compiler=mingw32 -Lc:\gsas\lib -lgsas
mv %1.pyd ..\binwin2.6
