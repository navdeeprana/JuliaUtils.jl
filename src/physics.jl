"""
Compute nearest neighbour distance for a given set of points in a 2D periodic box.
Useful for finding topological nearest neighbours.
$(TYPEDSIGNATURES)
"""
function nearest_neighbour_distance(x::AbstractArray{T}, y::AbstractArray{T}, L::T) where {T<:Real}
   N = length(x)
   d = zero(x)
   hL = 0.5 * L
   @inbounds for i in 1:N
      nnd = typemax(T)
      @inbounds for j in 1:N
         if (i == j)
               continue
         end
         dx = periodic_distance(abs(x[i] - x[j]), L, hL)
         dy = periodic_distance(abs(y[i] - y[j]), L, hL)
         dr = dx^2 + dy^2
         if (dr < nnd)
               nnd = dr
         end
      end
      d[i] = sqrt(nnd)
   end
   return d
end

"""
Compute radial distribution function for a given set of points in a 2D periodic box.
$(TYPEDSIGNATURES)
"""
function radial_distribution_function(x::Array{T}, y::Array{T}, L::T; nbin = 256) where {T<:Real}
   N = length(x)
   hL = 0.5 * L
   dbin = L / nbin
   r = dbin .* linspace(0, nbin, nbin)
   gr = zeros(nbin)
   @inbounds for i in 1:N
      @inbounds for j in i+1:N
         dx = periodic_distance(abs(x[i] - x[j]), L, hL)
         dy = periodic_distance(abs(y[i] - y[j]), L, hL)

         dr = sqrt(dx^2 + dy^2)
         ir = Int(ceil(dr / dbin))
         gr[ir] = gr[ir] + 1
      end
   end
   rho = (1 / (L^2)) * (N * (N - 1) / 2)
   @. gr = gr / (rho * 2 * π * r * dbin)
   return r, gr
end

"""
$(TYPEDSIGNATURES)
Compute order parameter given cartesian components of the field.\\
Args:\\

   - `u::Vararg{Array{T}}` [cartesian components of the field, accepts arrays of `Float64` or `ComplexF64`]\\

Examples:\\

 - for 3D vector field ` 𝐩 = (u,v,w)`; call the function as `order_parameter(u,v,w)` \\
 - for tensorial field `𝐐 = [u v; w z]`;  call the function as `order_parameter(u,v,w,z)` \\
 - Output is a tuple of `(|𝐦|, [⟨u⟩,⟨v⟩,⟨w⟩...])`, where `|𝐦|` is Frobenius norm of input tensor `𝐐` or vector `𝐩`.\\
"""
function order_parameter(u::Vararg{Array{T}}) where {T<:Union{Float64,ComplexF64}}
   op = zeros(T,length(u))
   for (i, ui) in enumerate(u)
      op[i] = sum(ui)/length(ui)
   end
   return hypot(op...),op
end
# https://math.stackexchange.com/questions/1958882/on-the-generalization-of-polar-coordinates-for-n-dimensions