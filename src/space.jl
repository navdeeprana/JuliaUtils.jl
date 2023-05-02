"""
---
numpy like linspace.
$(TYPEDSIGNATURES)
Args:\\
- [Req] x0::Real    (start point)\\
- [Req] xn::Real    (end point)\\
- [Req] n::Integer  (# of points in stop-start range.)\\
Note:\\
^ endpoint false.\\
^ useful for generating lattice of points in periodic box where the right hand edge is the box imaginary boundary.
"""
function linspace(x0::Real, xn::Real, n::Integer)
   return LinRange(x0, xn, n + 1)[1:end-1]
end
"""
$(SIGNATURES)
"""
function logspace(x0::Real, xn::Real, n::Integer;base=10)
   return base.^(range(log(base,x0),stop=log(base,xn),length=n+1)[1:end-1])
end
#####################

"""
Mesh in cartesian coordinates by explictly passing x (row) and y (col) range vectors.
Ouput is produced in classic (x,y) format with (id = row+nx*col) rule from left bottom corner
$(TYPEDSIGNATURES)
Args:\\
[Req] x::AbstractVector{T} (x coordinates range, row)\\
[Req] y::AbstractVector{T} (y coordinates range, col. Use y=x if square box)\\
Examples:\\
cartesian_mesh(0:1:10,0:1:10)
"""
function cartesian_mesh(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
   return (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
end
"""
Mesh in cartesian coordinates
$(TYPEDSIGNATURES)
Args:\\
[Req] nx::Integer  (# of points in x axis)\\
[Req] lx::Real     (length of the box in x axis)\\
[Opt] ny = nx      (# of points in y axis)\\
[Opt] ly = lx      (length of the box in y axis)\\
[Opt] x0 = 0       (ofset coordinates in x axis)\\
[Opt] y0 = 0       (ofset coordinates in y axis)\\
Note:\\
@ endpoint false. coordinates ranges x=[-lx/2, lx/2-1], y=[-ly/2, ly/2-1]\\
@ useful for generating lattice of points in periodic box where the right hand edge is the box boundary.\\
@ output is tuple follwoing (id = row+nx*col) rule from left bottom corner
"""
function cartesian_mesh(nx::Integer, lx::Real; ny = nx, ly = lx, x0 = 0, y0 = 0)
   x = linspace(-lx / 2, lx / 2, nx) .- x0
   y = (linspace(-ly / 2, ly / 2, ny) .- y0)
      return  x,y'
end
#####################

"""
``Mesh in polar coordinates.``
$(TYPEDSIGNATURES)
``Args:``\\
[Req] nr::Integer (# of points in radial direction)
[Req] lr::Real;   (length of the radius)
[Opt] nθ = nr     (# of points in θ direction)
[Opt] lθ = 2π     (length of θ direction)
[Opt] r0 = 0      (ofset in radial direction; not implemented yet)
[Opt] θ0 = 0      (ofset in theta direction; not implemented yet)
"""
function polar_mesh(nr::Integer, lr::Real; nθ = nr, lθ = 2π, r0 = 0, θ0 = 0)
   r = linspace(r0, lr, nr)
   θ = linspace(θ0, lθ, nθ)
      return  r, θ'
end

"""
``Convert cartesian mesh to polar mesh.``
$(TYPEDSIGNATURES)
``Args:``\\
[Req] x::AbstractVector{T} (explicit x coordinates in classic form, row)\\
[Req] y::AbstractVector{T} (explicit y coordinates in classic form, col)\\
[Opt] classic::Bool=false\\
``Examples:``\\
cartesian_mesh(0:1:10,0:1:10)\\
"""
function polar_mesh(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real}
   r = hypot.(x, y)
   θ = (atan.(y, x))
   return  r, θ
end
#####################

"""
``Picks logspaced elements from given 1D array``
$(TYPEDSIGNATURES)
``Args:``\\
[Req] x::AbstractVector;   (input vector of ordered data)\\
[Opt] base = 1.2  (choose base )
"""
function logspacedPick(x::AbstractVector{T}; base = 1.2) where {T<:Real}
   N = log(length(x)) / log(base)
   n = unique(floor.(Int, base .^ (1:N)))
   return x[n],n
end
function logspacedPick!(x::AbstractVector{T}; base = 1.2) where {T<:Real}
   sort!(x) #mutating orginal input vector
   N = log(length(x)) / log(base)
   n = unique(floor.(Int, base .^ (1:N)))
   return x[n],n
end
"""
``Picks logspaced radial elements from given 2D array after sorting the radius data.``\\
``Output is in radius^2 format``
$(TYPEDSIGNATURES)
``Args:``\\
[Req] x::AbstractVector;   (input vector x cords; explicit form)\\
[Req] y::AbstractVector;   (input vector y cords; explicit form)\\
[Opt] base = 1.2  (choose base)\\
"""
function logspacedPick(x::AbstractVector{T},y::AbstractVector{T}; base = 1.2) where {T<:Real}
   # r = @. (x^2 + y^2)^(0.5)
   r = @. (x^2 + y^2)
   sort!(r)
   return logspacedPick(r,base=base)
end
