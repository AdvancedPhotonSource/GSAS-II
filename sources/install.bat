set OUTDIR=%USERPROFILE%\.GSASII\bin
copy /y "%MESON_BUILD_ROOT%\sources\"*.so ${OUTDIR}\
copy /y "%MESON_BUILD_ROOT%\sources\"*/*.so ${OUTDIR}\
copy /y "%MESON_BUILD_ROOT%\sources\"LATTIC ${OUTDIR}\
copy /y "%MESON_BUILD_ROOT%\sources\"convcell ${OUTDIR}\
copy /y "%MESON_BUILD_ROOT%\sources\"GSAS*.txt ${OUTDIR}\
