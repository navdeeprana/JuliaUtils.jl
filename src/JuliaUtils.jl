"""
Main module for `JuliaUtils.jl` -- my julia language utility functions.

# Exports
$(EXPORTS)
"""
module JuliaUtils

export
    read_data,
    read_xydata,
    write_data,
    write_xydata

include("io.jl")

end
