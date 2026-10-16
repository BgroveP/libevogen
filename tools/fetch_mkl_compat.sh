#!/usr/bin/env bash
set -euo pipefail
d="/mnt/c/Users/au488376/OneDrive - Aarhus universitet/Documents/src/libevogen_jll/tools/mkl_compat"
b=https://raw.githubusercontent.com/Reference-LAPACK/lapack/master
tmp=$(mktemp -d)
cd "$tmp"
apt-get download liblapacke-dev >/dev/null 2>&1
dpkg -x liblapacke-dev*.deb pkg
find pkg -name "*.h" | xargs ls -la
cp pkg/usr/include/lapacke.h pkg/usr/include/lapacke_64.h \
   pkg/usr/include/lapacke_config.h pkg/usr/include/lapacke_mangling.h "$d"/
curl -fsSL -o "$d/cblas.h" "$b/CBLAS/include/cblas.h"
cat > "$d/mkl_lapacke.h" <<'EOF'
/* Minimal stand-in for Intel MKL's mkl_lapacke.h for building against
 * libevogen's headers without a oneAPI installation. The LAPACKE C
 * interface uses the same LAPACKE_* function names in MKL and in the
 * reference implementation, so the netlib lapacke.h provides identical
 * declarations. At runtime the symbols resolve against MKL (via
 * MKL_jll's libmkl_rt.so, see the JLL wrapper).
 *
 * Build-time header only; not needed by consumers of the artifact. */
#ifndef MKL_LAPACKE_COMPAT_H
#define MKL_LAPACKE_COMPAT_H
#include "lapacke.h"
#endif
EOF
ls -la "$d"
