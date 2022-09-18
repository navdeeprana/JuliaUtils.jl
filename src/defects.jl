"""
Generate a topological defect vector field in two dimensions.

$(SIGNATURES)
"""
function defect_vector_field(m::Integer, lx, nx; θ0 = 0, kwargs...)
    x, y = cartesian_mesh(nx, lx; kwargs...)
    θ = @. atan(y, x)
    u = @. cos(m * θ + θ0)
    v = @. sin(m * θ + θ0)
    return u, v
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
Count all topological zeros for a 2D orientation field.

$(SIGNATURES)
"""
function count_topological_charges(θ::Array{T,2}; threshold = 0.8 * (2π)) where {T<:Real}
    nx, ny = size(θ)
    x, y, c = Int[], Int[], Int[]

    @inbounds for j = 1:ny
        jnext = _get_next(j, ny)
        @inbounds for i = 1:nx
            inext = _get_next(i, nx)
            t1 = θ[i, j]
            t2 = θ[inext, j]
            t3 = θ[inext, jnext]
            t4 = θ[i, jnext]

            charge = _area(t2, t1) + _area(t3, t2) + _area(t4, t3) + _area(t1, t4)

            if abs(charge) > threshold
                push!(x, i)
                push!(y, j)
                push!(c, round(Int, charge / (2π)))
            end
        end
    end
    return x, y, c
end

function count_topological_charges(u::Array{T,2}, v::Array{T,2}; kwargs...) where {T<:Real}
    θ = @. atan(v, u)
    return count_topological_charges(θ; kwargs...)
end

function count_topological_charges(ϕ::Array{Complex{T},2}; kwargs...) where {T<:Real}
    θ = @. angle(ϕ)
    return count_topological_charges(θ; kwargs...)
end
