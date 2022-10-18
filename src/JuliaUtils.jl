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
    count_topological_charges,
    nearest_neighbour_distance,
    radial_distribution_function,
    @repeat

using DocStringExtensions

include("io.jl")
include("space.jl")
include("defects.jl")
include("physics.jl")
include("macros.jl")

end

module PyPlotUtils

export colors,
       figax,
       noticks!,
       set_fontsize!,
       set_formatter!

include("pyplot.jl")

end
