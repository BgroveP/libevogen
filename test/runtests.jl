using Test
using libevogen_jll
using Libdl

@testset "libevogen_jll" begin
    @test libevogen_jll.is_available() == (Sys.islinux() && Sys.ARCH === :x86_64)

    if !libevogen_jll.is_available()
        @warn "Platform not supported by libevogen_jll; skipping artifact checks" platform = (Sys.KERNEL, Sys.ARCH)
    else
        # exported products: paths usable in ccall / dlopen
        @test isfile(libevogen_jll.libevogen)
        @test libevogen_jll.libevogen == libevogen_jll.get_libevogen_path()
        @test isfile(libevogen_jll.get_libplinkio_path())
        @test isfile(libevogen_jll.libevogen_wrapper)
        @test libevogen_jll.libevogen_wrapper == libevogen_jll.get_libevogen_wrapper_path()

        # libevogen and libplinkio were dlopen'ed by __init__
        @test libevogen_jll.libevogen_handle != C_NULL
        @test libevogen_jll.libplinkio_handle != C_NULL

        # headers shipped in the artifact (needed to build wrappers)
        @test isdir(libevogen_jll.include_dir())
        @test isfile(joinpath(libevogen_jll.include_dir(), "Population.hpp"))

        # all libraries shipped, in the artifact lib dir
        libs = ["libintlc.so.5", "libimf.so", "libsvml.so", "libirng.so",
                "libiomp5.so", "libplinkio.so", "libevogen.so",
                "libevogen_wrapper.so"]
        for lib in libs
            @test isfile(joinpath(libevogen_jll.lib_dir(), lib))
        end

        # a real C++ symbol resolves from the loaded handle
        @test Libdl.dlsym(libevogen_jll.libevogen_handle,
                          :_ZN6evogen5TraitC2Ev) != C_NULL

        # the libcxxwrap wrapper exposes the module entry point. It is not
        # preloaded by the JLL (its libcxxwrap symbols resolve only once
        # the consumer has loaded libcxxwrap_julia), so dlopen it here
        # with RTLD_LAZY: function references stay unresolved until first
        # call, which is fine for this check.
        wh = Libdl.dlopen(libevogen_jll.libevogen_wrapper,
                          Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
        @test wh != C_NULL
        @test Libdl.dlsym(wh, :define_julia_module) != C_NULL
    end
end
