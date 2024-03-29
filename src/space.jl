"""
numpy like linspace.

$(SIGNATURES)
"""
linspace(x0::Real, xn::Real, n::Integer) = LinRange(x0, xn, n + 1)[1:end-1]

"""
Mesh in cartesian coordinates.

$(SIGNATURES)
"""
function cartesian_mesh(nx::Integer, lx::Real; ny = nx, ly = lx, x0 = 0, y0 = 0)
    x = linspace(-lx / 2, lx / 2, nx) .- x0
    y = linspace(-ly / 2, ly / 2, ny)' .- y0
    return x, y
end

"""
Mesh in polar coordinates.

$(SIGNATURES)
"""
function polar_mesh(nr::Integer, lr::Real; nθ = nr, lθ = 2π, r0 = 0, θ0 = 0)
    r = linspace(r0, lr, nr)
    θ = linspace(θ0, lθ, nθ)'
    return r, θ
end

"""
Convert cartesian mesh to polar mesh.

$(SIGNATURES)
"""
function polar_mesh(x::AbstractVector{T}, y::AbstractMatrix{T}) where {T<:Real}
    r = hypot.(x, y)
    θ = atan.(y, x)
    return r, θ
end

@inline periodic_distance(x::T, l::T, hl::T) where {T<:Real} = @. ifelse(abs(x) > hl, l - x, x)

@inline periodic_distance(x::T, l::T) where {T<:Real} = periodic_distance(x, l, 0.5 * l)

"""
Convert periodic data to an aperiodic one.

$(SIGNATURES)
"""
function aperiodic(x::Vector{T}, p = 2π) where {T<:Real}
    y = zero(x)
    x0, Δ, Δp = 0, 0, 0.9 * p
    for (i, xi) in enumerate(x)
        Δx = xi - x0
        if (abs(Δx) >= Δp)
            Δ = Δ - sign(Δx) * p
        end
        y[i] = xi + Δ
        x0 = xi
    end
    return y
end

"""
Pick logspaced elements from an array.

$(SIGNATURES)
"""
function logspaced(x::AbstractVector; base = 1.2)
    N = log(length(x)) / log(base)
    n = unique(floor.(Int, base .^ (1:N)))
    return x[n]
end
