path=%path:MinGW32=MinGW64%
if [%2] == [] goto ONE
	cd %2
	gfortran -c *.for -w -O2 -fno-automatic -finit-local-zero -malign-double -mwindows
        ar crs  lib%2.a *.o
	del *.o
	move  /Y lib%2.a ..
	cd ..
	f2py -c -m %1 %1.for --compiler=mingw32 -L./ -l%2
	goto TWO
:ONE
f2py -c -m %1 %1.for --compiler=mingw32 
:TWO
path=%path:MinGW64=MinGW32%
move %1.pyd ..\binwin64-2.6

