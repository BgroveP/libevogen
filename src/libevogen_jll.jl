# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libevogen_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libevogen")
JLLWrappers.@generate_main_file("libevogen", Base.UUID("923945fe-9270-4730-950a-f5510da5f31c"))
end  # module libevogen_jll
