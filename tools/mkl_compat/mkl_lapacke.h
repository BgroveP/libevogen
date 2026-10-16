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
