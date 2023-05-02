using Interpolations
"""
Interpolate a 2D data 
$(TYPEDSIGNATURES)
- [Req] data::AbstractArray{T,1} (input data vector, row ordered (0,0) at bottom left corner)
- [Req] order::T=1.3   (order of interpolation)
Note:\\
data is strutred with (row+nx*col) rule starting from bottom left corner.
"""
function Interp2D(data::AbstractArray{T,1}; order::T=1.3) where T<:Real
    
   IC = CubicSplineInterpolation((axes(data,1), axes(data,2)), data)

   finerx = LinRange(firstindex(data,1), lastindex(data,1), size(data,1) * order)
   finery = LinRange(firstindex(data,2), lastindex(data,2), size(data,2) * order)
   nx = length(finerx)
   ny = length(finery)

   data_interp = Array{Float64}(undef,nx,ny)
   for i ∈ 1:nx, j ∈ 1:ny
      data_interp[i,j] = IC(finerx[i],finery[j])
   end

   return finerx, finery, data_interp

end
