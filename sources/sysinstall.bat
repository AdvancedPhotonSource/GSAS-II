set OUTDIR=%MESON_SOURCE_ROOT%\GSASII\bin
mkdir %OUTDIR%
copy /y "%MESON_BUILD_ROOT%\sources\"*.pyd %OUTDIR%\
copy /y "%MESON_BUILD_ROOT%\sources\k_vec_cython\"*.pyd %OUTDIR%\
copy /y "%MESON_BUILD_ROOT%\sources\"LATTIC.exe %OUTDIR%\
copy /y "%MESON_BUILD_ROOT%\sources\"convcell.exe %OUTDIR%\
copy /y "%MESON_BUILD_ROOT%\sources\"GSAS*.txt %OUTDIR%\
copy /y %CONDA_PREFIX%\Library\bin\libgfortran.dll %OUTDIR%\
copy /y %CONDA_PREFIX%\Library\bin\libquadmath*.dll %OUTDIR%\
copy /y %CONDA_PREFIX%\Library\bin\libgcc*.dll %OUTDIR%\
REM copy /y %CONDA_PREFIX%\Library\bin\libwinthread*.dll %OUTDIR%\
