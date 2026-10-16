# Build the libcxxwrap wrapper (libevogen_wrapper.so) into the staged
# artifact (artifact_staging/lib/), against:
#   - the evogen C++ headers staged in artifact_staging/include,
#   - the libcxxwrap-julia headers from libcxxwrap_julia_jll,
#   - the Julia headers of the running Julia.
#
# Usage (run with a Julia whose environment has libcxxwrap_julia_jll
# and MKL_jll, at the same versions the consumers resolve):
#   julia --project=<env> tools/build_wrapper.jl [jlinterface.cpp]
#
#   jlinterface.cpp  path of the wrapper source
#                   (default: src/JLinterface.cpp in this repository)
#
# Run after tools/rebuild_libs.sh and before tools/install_artifact.jl.
#
# The wrapper is ABI-coupled to the libcxxwrap version and (via the
# jlcxx headers) to the Julia version family: rebuild it and re-install
# the artifact whenever libcxxwrap_julia_jll gets an ABI-breaking bump,
# and keep the compat entries of the consuming package (Evogen.jl) in
# sync with what was used here.
using libcxxwrap_julia_jll
using MKL_jll

const PKG_ROOT = dirname(@__DIR__)
const STAGING = joinpath(PKG_ROOT, "artifact_staging")
const OUT = joinpath(STAGING, "lib", "libevogen_wrapper.so")
const CPP = length(ARGS) >= 1 ? ARGS[1] :
    normpath(joinpath(PKG_ROOT, "src", "JLinterface.cpp"))

isdir(joinpath(STAGING, "lib")) || error(
    "no staged artifact found; run tools/rebuild_libs.sh first")
isfile(CPP) || error("wrapper source not found: $CPP")

julia_inc = normpath(joinpath(dirname(Sys.BINDIR), "include", "julia"))
evogen_inc = joinpath(STAGING, "include")
jlcxx_inc = joinpath(libcxxwrap_julia_jll.artifact_dir, "include", "jlcxx")
evogen_lib = joinpath(STAGING, "lib")

# The evogen headers include <mkl_lapacke.h> and <cblas.h> but are built
# without a oneAPI installation; tools/mkl_compat/ provides the reference
# LAPACKE/CBLAS declarations (identical LAPACKE_*/cblas_* names to MKL's
# interface). libmkl_rt is linked with --as-needed, so it ends up as a
# NEEDED entry only if the bound methods reference cblas/LAPACKE symbols;
# if a rebuild does record that NEEDED entry, add MKL_jll as a dependency
# of the JLL and preload libmkl_rt in __init__ (libevogen.so itself has
# MKL statically linked and needs no MKL runtime).
mkl_compat_inc = joinpath(PKG_ROOT, "tools", "mkl_compat")
mkl_lib = joinpath(MKL_jll.artifact_dir, "lib")

@info "Compiling wrapper" cpp = CPP out = OUT
@info "Using libcxxwrap headers from" libcxxwrap_julia_jll.artifact_dir

cmd = `g++ -std=c++20 -fPIC
       -I$julia_inc -I$evogen_inc -I$jlcxx_inc -I$mkl_compat_inc
       -shared $CPP
       -L$evogen_lib
       -Wl,-rpath,\$ORIGIN
       -Wl,--no-as-needed
       -levogen -lplinkio
       -Wl,--as-needed
       -L$mkl_lib -lmkl_rt
       -o $OUT`

run(cmd)
@info "Wrapper built" out = OUT
println("Next step: julia tools/install_artifact.jl")
