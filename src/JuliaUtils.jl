"""
Main module for `JuliaUtils.jl` -- my julia language utility functions.

# Exports
$(EXPORTS)
"""
module JuliaUtils

using DocStringExtensions

export read_data, read_xydata, write_data, write_xydata

include("io.jl")

export linspace, logspaced, cartesian_mesh, polar_mesh, periodic_distance, aperiodic

include("space.jl")

export defect_vector_field, find_defects, charge

include("defects.jl")

export nearest_neighbour_distance, radial_distribution_function, @repeat

include("physics.jl")

export @repeat

include("macros.jl")

module MakieUtils

export figax

include("makie.jl")

end

using PrecompileTools

@setup_workload begin
    @compile_workload begin
        cartesian_mesh(8, 1.0)
        polar_mesh(8, 1.0)
    end
end

end
