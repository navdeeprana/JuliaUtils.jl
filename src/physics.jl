"""
Compute nearest neighbour distance for a given set of points in a 2D periodic box.

$(SIGNATURES)
"""
function nearest_neighbour_distance(x::Array{T}, y::Array{T}, L::T) where {T<:Real}
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

$(SIGNATURES)
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
