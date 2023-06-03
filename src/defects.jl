"""
Generate a single topological defect vector field in two dimensions.
$(TYPEDSIGNATURES) 
Args:\\
[Req] q::Integer   (value of the charge)\\
[Req] lx::Real     (linear size of the box. The square box is used.)\\
[Req] m::Integer   (value of topological charge.)\\
[Req] nx::Integer  (number of grid points in the box, with ny=nx)\\
[Opt] θ0 = 0       (constant overall global phase, default is zero.)\\
[Opt] kwargs...    (extra keyword arguments for cartesian_mesh() function.)
Note:\\
functionn generates its own box and cartesian mesh. Use other method to provide mesh coordinates explictly.
@Examples
a = 2
b =4 
a+b
"""
function defect_vector_field(q::Integer, lx::Real, nx::Integer; θ0 = 0, kwargs...)
   x, y = cartesian_mesh(nx, lx; kwargs...)
   θ = @. atan(y, x)
   u = @. cos(q * θ + θ0)
   v = @. sin(q * θ + θ0)
   return u, v
end
"""
Generate vector field due to many topological defectsin two dimensions.
$(TYPEDSIGNATURES) 
Args:\\
[Req] x::AbstractArray{U} (x coordinates, in XY Form, xy' Form or mesh Form)\\
[Req] y::AbstractArray{U} (x coordinates, in XY Form, xy' Form or mesh Form)\\
[Req] Q::AbstractVector{<:Real} (Array of charge values,[q1,q2,...])\\
[Req] postQ::AbstractVector{T}  (Array of tuples for coordinates of charges, [(ax1,ay1),(ax2,ay2),...])\\
[Opt] θ0 = 0 (global phase for the vector field)
"""
function defect_vector_field(
   x::AbstractArray{U},
   y::AbstractArray{U},
   Q::AbstractVector{<:Real},
   postQ::AbstractVector{T};
   θ0 = 0) where {T<:Tuple{Real,Real},U<:Real}

   @assert (typeof(x)<:AbstractVector && y'==adjoint(x) ? false : true)  "expected (x,y <: AbstractVector{Real} in XY form) or (x,y <:AbstractArray{Real,2} in mesh form) or (y<:Adjoint of x)"
   @assert (size(Q)[1] == size(postQ)[1]) "Number of charges does not have equivalent postion array."

   θ = fill!(similar(x.*y),0.0)
   @inbounds for (q,(ax,ay)) in zip(Q,postQ)
     @. θ +=  q*atan(y-ay, x-ax)
   end
   u = @. cos(θ + θ0)
   v = @. sin(θ + θ0)
   return u,v
end

@inline function _get_next(i::T, N::T) where {T<:Integer}
   i == N ? 1 : i + 1
end

@inline function _area(tn::T, t0::T) where {T<:Real}
   dt = tn - t0
   if (dt > π)
      dt = dt - 2π
   end
   if (dt < -π)
      dt = dt + 2π
   end
   return dt
end

"""
$(METHODLIST)
"""
function find_defects end
"""
Count all topological zeros for a 2D orientation field.
$(TYPEDSIGNATURES)
Args:\\
[Req] θ::AbstractArray{T,2}; (Array of θ values, make sure they are column ordered)\\
[Opt] θ0::Real = 0.8 * (2π) (optional threshold value above which a charge is counted. Doesn't affect much.)\\
[Opt] mx::Integer = 0 (scan in a small window of the box, anchored at (0,0) and streched upto (nx-mx,ny-my))\\
[Opt] my::Integer = mx (scan in a small window of the box, anchored at (0,0) and streched upto (nx-mx,ny-my))
"""
function find_defects(θ::AbstractArray{T,2}; θ0::Real = 0.8 * (2π), mx::Integer = 0, my::Integer = mx) where {T<:Real}
   nx, ny = size(θ)
   x, y, c = Int[], Int[], Int[]

   @inbounds for j in 1:ny-my
      jnext = _get_next(j, ny)
      @inbounds for i in 1:nx-mx
         inext = _get_next(i, nx)
         t1 = θ[i, j]
         t2 = θ[inext, j]
         t3 = θ[inext, jnext]
         t4 = θ[i, jnext]

         Δθ = _area(t2, t1) + _area(t3, t2) + _area(t4, t3) + _area(t1, t4)

         if abs(Δθ) > θ0
               push!(x, i)
               push!(y, j)
               push!(c, round(Int, Δθ / (2π)))
         end
      end
   end
   return x, y, c
end
"""
Count all topological zeros for a 2D orientation field given cos and sin values.
$(TYPEDSIGNATURES)
Args:\\
[Req] u::AbstractArray{T,2}; (Array of cos(θ) values, make sure they are column ordered)\\
[Req] u::AbstractArray{T,2}; (Array of sin(θ) values, make sure they are column ordered)\\
[Opt] kwargs... (Same as find_defects() function)\\
"""
function find_defects(u::AbstractArray{T,2}, v::AbstractArray{T,2}; kwargs...) where {T<:Real}
   θ = @. atan(v, u)
   return find_defects(θ; kwargs...)
end
"""
Count all topological zeros for a 2D orientation field given complex field values.
$(TYPEDSIGNATURES)
Args:\\
[Req] ϕ::AbstractArray{Complex{T},2} (Array of cos(θ) + ι sin(θ), column ordered)\\
[Opt] kwargs... (Same as find_defects() function)\\
"""
function find_defects(ϕ::AbstractArray{Complex{T},2}; kwargs...) where {T<:Real}
   θ = @. angle(ϕ)
   return find_defects(θ; kwargs...)
end

"""
Compute topological charge at points (x,y) for the orientation field θ
$(TYPEDSIGNATURES)
"""
function chargeAt(x, y, θ::AbstractArray{T,2}; Δ = 1) where {T<:Real}
   nx, ny = size(θ)
   c = zeros(Float64, length(x))
   for (n, (xn, yn)) in enumerate(zip(x, y))
      for Δy in -Δ:+Δ
         for Δx in -Δ:+Δ
               i, j = mod(xn + Δx - 1, nx) + 1, mod(yn + Δy - 1, ny) + 1
               inext = _get_next(i, nx)
               jnext = _get_next(j, ny)
               t1 = θ[i, j]
               t2 = θ[inext, j]
               t3 = θ[inext, jnext]
               t4 = θ[i, jnext]
               c[n] += (_area(t2, t1) + _area(t3, t2) + _area(t4, t3) + _area(t1, t4)) / (2π)
         end
      end
   end
   return c
end
"""
Compute topological charge at points (x,y) for the vector field cos(θ),sin(θ)
$(TYPEDSIGNATURES)
"""
function chargeAt(x, y, u::AbstractArray{T,2}, v::AbstractArray{T,2}; kwargs...) where {T<:Real}
   θ = @. atan(v, u)
   return chargeAt(x, y, θ; kwargs...)
end
"""
Compute topological charge at points (x,y) for the complex field cos(θ)+ ι sin(θ)
$(TYPEDSIGNATURES)
"""
function chargeAt(x, y, ϕ::AbstractArray{Complex{T},2}; kwargs...) where {T<:Real}
   θ = @. angle(ϕ)
   return chargeAt(x, y, θ; kwargs...)
end