"""
Wraps angle θ in [-π,π) or [0,2π) domain.\\
$(TYPEDSIGNATURES)
[Req] θ::Real (input angle value)
[Opt] domain=π (domain type, π => [-π,π) and 2π =  [0,2π) domain)
"""
function wrap_angle(θ::Real;domain=π)
   @assert domain==π || domain==2π "'domain' can take values of π or 2π, default is π"
   θ = mod(θ + domain,2π)
   if (θ < 0)
     θ += 2π
   end
   return θ - domain
end
"""
Wraps angle difference (θ2 - θ1) in [-π,π) or [0,2π) domain.\\
benchmark with _area(θ2,θ1) before using this function.
$(TYPEDSIGNATURES)
[Req] θ1::Real (input angle1 value)
[Req] θ2::Real (input angle2 value)
[Opt] domain=π (domain type, π => [-π,π) and 2π =  [0,2π) domain)
"""
function wrap_angle(θ2::T,θ1::T;domain=π) where {T<:Real}
   diff = θ2-θ1
   return anglewrap(diff,domain=domain)
end
# function anglewrap(angle::Real)
#    angle = mod(angle + π,2π);
#    if (angle < 0)
#       angle += 2π;
#    end
#    return angle - π;
# end

#####################
"""
find periodic distance provided absolute distance abs(xi - xj).
$(TYPEDSIGNATURES)
Args:\\
[Req] x::T,  (absolute distance)
[Req] l::T,  (linear box size)
[Req] hl::T  (half length of box size)
"""
@inline function periodic_distance(x::Real, l::Real, hl::Real)
   @assert !(x < 0.0) "provide absolute distance between two coordinates abs(xj - xi)"
   @. ifelse(abs(x) > hl, l - x, x)
end
"""
find periodic distance provided absolute distance abs(xi - xj).
$(TYPEDSIGNATURES)
Args:\\
[Req] x::T,  (absolute distance)\\
[Req] l::T,  (linear box size)\\
"""
@inline function periodic_distance(x::Real, l::Real)
   @assert !(x<0.0) "provide absolute distance between two coordinates abs(xj - xi)"
   return periodic_distance(x, l, 0.5 * l)
end
"""
Wrap particle coordinate in [-hl,hl) open range with origin at center.
$(TYPEDSIGNATURES)
Args:\\
[Req] x::T,  (absolute distance)\\
[Req] l::T,  (linear box size)\\
[Req] hl::T  (half length of box size)
"""
@inline function wrap_cords(x::Real, l::Real, hl::Real)
   if (abs(x)>hl)
      if (x >= hl)
         x = x - l
      end
      if (x < -hl)
         x = x + l
      end
   end
   return x
end
"""
Wrap particle coordinate in [-hl,hl) open range with origin at center.
$(TYPEDSIGNATURES)
Args:\\
[Req] x::T,  (absolute distance)\\
[Req] l::T,  (linear box size)
"""
@inline function wrap_cords(x::Real, l::Real)
   return wrap_cords(x,l,0.5*l)
end
# """
# Wrap particle coordinate in [-hl,hl] closed range with origin at center.
# $(TYPEDSIGNATURES)
# Args:\\
# [Req] x::T,  (absolute distance)\\
# [Req] l::T,  (linear box size)
# """
# function periodic_cord(x::T, l::T, hl::T) where {T<:Real}
#    if abs(x) > hl
#       x -= l * round(x/l)
#    end
#    return x
# end
#############################################################
"""
Convert periodic data to an aperiodic one.
$(TYPEDSIGNATURES)
Args:\\
[Req] x::Vector{T};  (input data vector)\\
[Opt] period::Real = 2π
"""
function aperiodic(x::AbstractVector{T}; period::Real = 2π) where {T<:Real}
   y = zero(x)
   x0, Δ, Δp = 0, 0, 0.9 * period
   for (i, xi) in enumerate(x)
      Δx = xi - x0
      if (abs(Δx) >= Δp)
         Δ = Δ - sign(Δx) * period
      end
      y[i] = xi + Δ
      x0 = xi
   end
   return y
end
#####################


