# libevogen_jll

Julia JLL package wrapping the prebuilt `libevogen` shared library, its
runtime dependencies, the evogen C++ headers, and a prebuilt libcxxwrap
wrapper (`libevogen_wrapper.so`) that Julia packages can load directly
with [CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl).

The package follows the standard JLL anatomy
([docs.binarybuilder.org/stable/jll/](https://docs.binarybuilder.org/stable/jll/)):
an `Artifacts.toml` binding a `libevogen` artifact (linux x86_64, glibc),
a `baremodule` main file generated via `JLLWrappers` macros, and a
platform wrapper in `src/wrappers/`. The only supported platform is
`x86_64-linux-gnu`; on other platforms the package loads fine,
`is_available()` returns `false` and the library variables stay unset.

```
libevogen_jll/
├── Project.toml            # deps: JLLWrappers, Libdl, Artifacts, Preferences
├── Artifacts.toml          # [[libevogen]] binding (git-tree-sha1, no download)
├── artifacts/
│   └── libevogen/          # the artifact itself, committed in-repo
│       ├── lib/            # .so files (incl. libevogen_wrapper.so)
│       └── include/        # evogen C++ headers
├── lib/, include/          # pristine upstream binaries/headers (rebuild inputs)
├── .gitignore
├── src/
│   ├── libevogen_jll.jl           # main file: platform dispatch
│   ├── JLinterface.cpp            # wrapper source (compiled into the artifact)
│   └── wrappers/
│       └── x86_64-linux-gnu.jl     # products + library loading
├── test/runtests.jl
└── tools/                  # maintainer build scripts, not needed at runtime
    ├── rebuild_libs.sh            # stage patched lib/ + include/
    ├── build_wrapper.jl           # compile libevogen_wrapper.so
    ├── install_artifact.jl        # publish in-repo artifact + bind + install
    ├── fetch_mkl_compat.sh        # refresh tools/mkl_compat/ build headers
    └── mkl_compat/                # vendored netlib LAPACKE/CBLAS headers
```

## Exports and API

Standard JLL surface:

- `libevogen`: absolute path of `libevogen.so` (usable in `ccall`).
- `libplinkio`: absolute path of `libplinkio.so`.
- `libevogen_wrapper`: absolute path of `libevogen_wrapper.so`, the
  prebuilt libcxxwrap wrapper to pass to `CxxWrap.@wrapmodule`.
- `is_available()`: `true` iff the current platform is supported.
- Unexported, as in every JLL: `artifact_dir`, `PATH`, `LIBPATH`,
  `get_libevogen_path()`, `libevogen_handle`, `libplinkio_handle`,
  `get_libevogen_wrapper_path()`, `dev_jll()`.

Additions for wrapper-building packages:

- `lib_dir()`: the artifact's `lib/` directory.
- `include_dir()`: the artifact's `include/` directory (evogen C++ headers).

## The artifact

The artifact content is committed in this repository under
`artifacts/libevogen/` and bound in `Artifacts.toml` by its
`git-tree-sha1`. Pkg never installs such an artifact into the depot: it
silently skips bindings without a download section, so the platform
wrapper's `find_artifact_dir()` falls back to the copy committed in the
package's own git clone. No download URL and no public release needed,
so it works for private repositories too.

It contains:

- `lib/` - `libevogen.so` and its dependencies (the Intel oneAPI runtime
  libraries `libimf`, `libintlc.so.5`, `libiomp5`, `libirng`, `libsvml`,
  and `libplinkio`), plus `libevogen_wrapper.so`;
- `include/` - the evogen C++ headers.

`libevogen.so` and `libplinkio.so` are built without a SONAME and without
an rpath, so out of the box they only find their sibling libraries via
`LD_LIBRARY_PATH` or oneAPI env vars. The artifact ships them patched with
a `$ORIGIN` runpath and a matching SONAME (see `tools/rebuild_libs.sh`),
so they resolve their siblings wherever they are installed and register
under their names for NEEDED lookup. The JLL `__init__` additionally
preloads the Intel runtimes with `RTLD_GLOBAL`, so symbol resolution never
depends on the environment - no oneAPI `setvars.sh` needed.

MKL is *not* bundled: `libevogen.so` has MKL statically linked inside it.
The wrapper is compiled against the vendored reference LAPACKE/CBLAS
declarations in `tools/mkl_compat/` (build-time only) and linked against
`MKL_jll`'s `libmkl_rt.so` with `--as-needed`; the methods bound by the
current wrapper have no BLAS/LAPACK references, so the linker drops
`libmkl_rt` and consumers need nothing beyond this artifact. If future
wrapper bindings reach the cblas/LAPACKE code paths, add `MKL_jll` as a
dependency and preload `libmkl_rt` in `__init__` (see the comment in
`tools/build_wrapper.jl`).

`libevogen_wrapper.so` is deliberately **not** preloaded by the JLL: it
has undefined references to `libcxxwrap_julia` symbols, which only resolve
once the consumer (CxxWrap) has loaded `libcxxwrap_julia`. The JLL exports
its path and the consumer dlopens it.

## Installing

The artifact lives in this repository, so installation is just adding
the package (works while the repository is private):

```julia
using Pkg
Pkg.add(url = "https://github.com/BgroveP/libevogen_jll")
# or, for development: Pkg.develop(path = "libevogen_jll")
```

No manual step is needed on any machine: the JLL resolves the binaries
from the `artifacts/libevogen/` tree committed in its own git clone, and
prefers the depot artifact when one was installed (e.g. by
`tools/install_artifact.jl`). The package must be installed through git
(or otherwise come with its `artifacts/libevogen/` directory), because
the fallback reads that directory.

## Using the wrapper from a Julia package (e.g. Evogen.jl)

```julia
module MyPackage

using CxxWrap          # loads libcxxwrap_julia first
using libevogen_jll    # loads the artifact; __init__ preloads libevogen

@wrapmodule(() -> libevogen_jll.libevogen_wrapper)

function __init__()
    @initcxx
end

end
```

`CxxWrap` must be loaded before `@wrapmodule` is expanded so that the
wrapper's `libcxxwrap_julia` symbols resolve.

### ABI pinning

`libevogen_wrapper.so` is built against the headers of one
`libcxxwrap_julia_jll` version (currently 0.14.x) and one Julia version
family (currently 1.13). Consuming packages must pin compatible versions
in their `[compat]`:

```toml
[compat]
CxxWrap = "0.17"
libcxxwrap_julia_jll = "0.14"
julia = "1.13"
```

Bumping these requires rebuilding the wrapper: see the next section.

## Rebuilding the artifact (maintainer)

When new binaries of `libevogen.so` arrive, or when the wrapper must be
rebuilt (libcxxwrap_julia_jll or Julia version bump), from the package
root, on Linux:

```sh
tools/rebuild_libs.sh [source_root]      # default: lib/ + include/ in this repo
julia --project=<env with libcxxwrap_julia_jll + MKL_jll> tools/build_wrapper.jl
julia tools/install_artifact.jl          # publishes artifacts/libevogen/ + binds it
git add Artifacts.toml && git commit
```

To pick up new upstream binaries, replace `lib/` and `include/` with the
new release first (or pass its root as `source_root`). The repository is
otherwise self-contained: the wrapper source `src/JLinterface.cpp` is
vendored here.

`rebuild_libs.sh` copies the pristine binaries and headers from `lib/`
and `include/` into `artifact_staging/` and embeds the `$ORIGIN` rpath
and SONAME into `libevogen.so` and `libplinkio.so` - with `patchelf` if
installed, otherwise with a stdlib-only Python ELF patcher embedded in
the script (no other bytes are touched). `build_wrapper.jl` compiles
`src/JLinterface.cpp` against the staged headers, the
libcxxwrap headers from `libcxxwrap_julia_jll`, the running Julia's
headers, and the vendored reference LAPACKE/CBLAS declarations in
`tools/mkl_compat/` (build-time only; see the MKL paragraph above).
`install_artifact.jl` copies the staging tree to
`artifacts/libevogen/`, `git add`s it, binds the staged git tree hash in
`Artifacts.toml`, and installs the artifact into the local depot. None of
these tools are needed at runtime.

## Testing

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

On x86_64 Linux: the artifact loads, all library products resolve, and a
C++ symbol as well as the wrapper's `define_julia_module` entry point can
be found via `dlsym`. On other platforms the artifact checks are skipped.

## Licensing

All rights to the original libraries belong to their creators; see
the license overview in the [libevogen repository](https://github.com/BgroveP/libevogen)
(`libevogen`: V. Milkevych; Intel runtime libraries: Intel EULA;
`libplinkio`: M. Frånberg, BSD). The wrapper source (now
`src/JLinterface.cpp` in this repository, originally
`libevogen/src/JLinterface.cpp`) is by Bjarke Poulsen. The vendored headers
in `tools/mkl_compat/` come from the netlib LAPACK/LAPACKE and CBLAS
reference implementations (BSD-3) and the Debian `liblapacke-dev`
package, and are used only to build the wrapper.
