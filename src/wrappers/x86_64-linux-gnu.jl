# Wrapper script for libevogen_jll for x86_64-linux-gnu.
#
# Follows the structure of BinaryBuilder-generated JLL wrappers
# (see https://docs.binarybuilder.org/stable/jll/). The Intel oneAPI
# runtime libraries are preloaded with RTLD_GLOBAL in __init__ (see the
# comment there for why this JLL needs it).
export libevogen, libplinkio, libevogen_wrapper, lib_dir, include_dir

# Equivalent of JLLWrappers.@generate_wrapper_header("libevogen") with a
# fallback. The Artifacts.toml entry for libevogen has no download
# section (the binaries are committed in this repository instead), so
# Pkg silently skips it and a fresh machine never gets the artifact in
# its depot. Fall back to the copy committed next to this source file;
# the depot artifact keeps precedence when it was installed (e.g. via
# tools/install_artifact.jl).
function find_artifact_dir()
    try
        return @artifact_str("libevogen")
    catch
        committed = normpath(joinpath(@__DIR__, "..", "..", "artifacts", "libevogen"))
        isdir(committed) || rethrow()
        return committed
    end
end
eager_mode() = nothing
if ccall(:jl_generating_output, Cint, ()) == 1
    Base.precompile(find_artifact_dir, ())
    Base.precompile(eager_mode, ())
end

JLLWrappers.@declare_library_product(libevogen, "libevogen.so")
JLLWrappers.@declare_library_product(libplinkio, "libplinkio.so")
JLLWrappers.@declare_library_product(libevogen_wrapper, "libevogen_wrapper.so")

# Conveniences for the Julia interface package: the evogen C++ headers
# shipped in the artifact are needed to build a libcxxwrap wrapper
# against libevogen.
lib_dir() = joinpath(artifact_dir, "lib")
include_dir() = joinpath(artifact_dir, "include")

function __init__()
    JLLWrappers.@generate_init_header()

    # libevogen.so and libplinkio.so ship without an rpath and without a
    # SONAME, so out of the box they only find their sibling libraries
    # via LD_LIBRARY_PATH or oneAPI env vars. The binaries in this
    # artifact are patched with a $ORIGIN rpath and a matching SONAME
    # (tools/rebuild_libs.sh), so they resolve their siblings in lib/
    # and register under their own names for NEEDED lookup. Preload the
    # Intel oneAPI runtime libraries anyway, with RTLD_GLOBAL, so symbol
    # resolution never depends on search paths.
    for lib in ("libintlc.so.5", "libimf.so", "libsvml.so", "libirng.so", "libiomp5.so")
        Libdl.dlopen(joinpath(artifact_dir, "lib", lib),
                     Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
    end

    JLLWrappers.@init_library_product(
        libplinkio,
        "lib/libplinkio.so",
        RTLD_LAZY | RTLD_GLOBAL,
    )
    JLLWrappers.@init_library_product(
        libevogen,
        "lib/libevogen.so",
        RTLD_LAZY | Libdl.RTLD_GLOBAL,
    )

    # The libcxxwrap wrapper is NOT dlopen'ed here: it has undefined
    # references to libcxxwrap_julia symbols, which only resolve once the
    # consumer (e.g. CxxWrap) has loaded libcxxwrap_julia with
    # RTLD_GLOBAL. Expose the path and let the consumer dlopen it.
    JLLWrappers.@init_library_product(
        libevogen_wrapper,
        "lib/libevogen_wrapper.so",
        nothing,
    )

    JLLWrappers.@generate_init_footer()
end  # __init__()
