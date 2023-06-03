"""
Main module for `Julia Useful functions`
# Exports
$(EXPORTS)
"""
module JuliaUtils
   #core exported functions
   export read_data,
      read_xydata,
      write_data,
      write_xydata,
      linspace,
      logspace,
      logspacedPick,
      logspacedPick!,
      cartesian_mesh,
      polar_mesh,
      periodic_distance,
      aperiodic,
      defect_vector_field,
      find_defects,
      chargeAt,
      nearest_neighbour_distance,
      radial_distribution_function,
      Interp2D,
      order_parameter,
      @repeat

   using DocStringExtensions
   include("io.jl")
   include("space.jl")
   include("defects.jl")
   include("physics.jl")
   include("macros.jl")
   include("interpolations.jl")
   include("periodic.jl")
   include("spectrum.jl")
   include("pyplot.jl")
   include("makie.jl")

   # include all files recursively from the `src` folder
   # for (root, dirs, files) in walkdir("src")
   #    for file in files
   #       # check if the file is a Julia file
   #       if endswith(file, ".jl")
   #          # construct the path to the file
   #          path = joinpath(root, file)
   #          # include the file in the module
   #          Base.include(path)
   #       end
   #    end
   # end




   # module PyPlotUtils
   # export colors,
   #    figax,
   #    noticks!,
   #    fontsize!,
   #    formatter!
   #    include("pyplot.jl")
   # end
   
   # module MakieUtils
   # using DocStringExtensions
   # export figax
   #    include("makie.jl")
   # end

end
