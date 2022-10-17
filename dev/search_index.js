var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = JuliaUtils","category":"page"},{"location":"#JuliaUtils","page":"Home","title":"JuliaUtils","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for JuliaUtils.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [JuliaUtils]","category":"page"},{"location":"#JuliaUtils.JuliaUtils","page":"Home","title":"JuliaUtils.JuliaUtils","text":"Main module for JuliaUtils.jl – my julia language utility functions.\n\nExports\n\naperiodic\ncartesian_mesh\ncount_topological_charges\ndefect_vector_field\nlinspace\nnearest_neighbour_distance\nperiodic_distance\npolar_mesh\nradial_distribution_function\nread_data\nread_xydata\n@repeat\nwrite_data\nwrite_xydata\n\n\n\n\n\n","category":"module"},{"location":"#JuliaUtils.aperiodic-Union{Tuple{Vector{T}}, Tuple{T}, Tuple{Vector{T}, Any}} where T<:Real","page":"Home","title":"JuliaUtils.aperiodic","text":"Convert periodic data to an aperiodic one.\n\naperiodic(x)\naperiodic(x, p)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.cartesian_mesh-Tuple{Integer, Real}","page":"Home","title":"JuliaUtils.cartesian_mesh","text":"Mesh in cartesian coordinates.\n\ncartesian_mesh(nx, lx; ny, ly, x0, y0)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.count_topological_charges-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:Real","page":"Home","title":"JuliaUtils.count_topological_charges","text":"Count all topological zeros for a 2D orientation field.\n\ncount_topological_charges(θ; θ0, mx, my)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.defect_vector_field-Tuple{Integer, Any, Any}","page":"Home","title":"JuliaUtils.defect_vector_field","text":"Generate a topological defect vector field in two dimensions.\n\ndefect_vector_field(m, lx, nx; θ0, kwargs...)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.linspace-Tuple{Real, Real, Integer}","page":"Home","title":"JuliaUtils.linspace","text":"numpy like linspace.\n\nlinspace(x0, xn, n)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.nearest_neighbour_distance-Union{Tuple{T}, Tuple{Array{T}, Array{T}, T}} where T<:Real","page":"Home","title":"JuliaUtils.nearest_neighbour_distance","text":"Compute nearest neighbour distance for a given set of points in a 2D periodic box.\n\nnearest_neighbour_distance(x, y, L)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.polar_mesh-Tuple{Integer, Real}","page":"Home","title":"JuliaUtils.polar_mesh","text":"Mesh in polar coordinates.\n\npolar_mesh(nr, lr; nθ, lθ, r0, θ0)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.polar_mesh-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractMatrix{T}}} where T<:Real","page":"Home","title":"JuliaUtils.polar_mesh","text":"Convert cartesian mesh to polar mesh.\n\npolar_mesh(x, y)\n\n\n\n\n\n\n","category":"method"},{"location":"#JuliaUtils.radial_distribution_function-Union{Tuple{T}, Tuple{Array{T}, Array{T}, T}} where T<:Real","page":"Home","title":"JuliaUtils.radial_distribution_function","text":"Compute radial distribution function for a given set of points in a 2D periodic box.\n\n\n\n\n\n","category":"method"}]
}
