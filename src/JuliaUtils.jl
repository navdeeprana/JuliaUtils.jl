"""
Main module for `JuliaUtils.jl` -- my julia language utility functions.

# Exports
$(EXPORTS)
"""
module JuliaUtils

export read_data,
    read_xydata,
    write_data,
    write_xydata,
    linspace,
    cartesian_mesh,
    polar_mesh,
    periodic_distance,
    aperiodic,
    defect_vector_field,
    count_topological_charges

using DocStringExtensions

include("io.jl")
include("space.jl")
include("defects.jl")

end
