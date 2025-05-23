#!/bin/bash
OUTDIR=${MESON_SOURCE_ROOT}/GSASII/bin
mkdir -p ${OUTDIR}
cp -v "${MESON_BUILD_ROOT}/sources/"*.so ${OUTDIR}/
cp -v "${MESON_BUILD_ROOT}/sources/"*/*.so ${OUTDIR}/
cp -v "${MESON_BUILD_ROOT}/sources/"LATTIC ${OUTDIR}/
cp -v "${MESON_BUILD_ROOT}/sources/"convcell ${OUTDIR}/
cp -v "${MESON_BUILD_ROOT}/sources/"GSAS*.txt ${OUTDIR}/
