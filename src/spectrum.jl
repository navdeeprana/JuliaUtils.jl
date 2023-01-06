using AbstractFFTs: fftfreq, rfftfreq

struct Spectrum
    D  :: Integer
    L  :: Float64
    Δk :: Float64
    k  :: Vector{Float64}
    sk :: Vector{Float64}
end

function Spectrum(N, D; L = 2π)
    Δk = 2π / L
    k  = @. Δk * LinRange(0, N, N + 1)
    sk = zeros(N + 1)
    return Spectrum(D, L, Δk, k, sk)
end

struct Correlator
    D  :: Integer
    r  :: Vector{Float64}
    cr :: Vector{Float64}
end

function Correlator(N, D; L = π)
    r  = LinRange(0, L, N + 1)
    cr = zeros(N + 1)
    return Correlator(D, r, cr)
end


@inline function _kmax(kiter)
    return round(Int, hypot(maximum(kiter)...)) + 1
end

function _kgrid(nx::T, xfftfreq::Function) where {T<:Integer}
    kiter = Iterators.product(xfftfreq(nx))
    return kiter, _kmax(kiter)
end

function _kgrid(nx::T, ny::T, xfftfreq::Function) where {T<:Integer}
    kiter = Iterators.product(xfftfreq(nx), fftfreq(ny, ny))
    return kiter, _kmax(kiter)
end

function _kgrid(nx::T, ny::T, nz::T, xfftfreq::Function) where {T<:Integer}
    kiter = Iterators.product(xfftfreq(nx), fftfreq(ny, ny), fftfreq(nz, nz))
    return kiter, _kmax(kiter)
end

@inline function _kbin(k::Tuple{Vararg{T}}) where {T<:Real}
    return round(Int, hypot(k...)) + 1
end

function spectrum!(S::Spectrum, kiter::Iterators.ProductIterator, uk::Array{ComplexF64}; preserve = true)
    N = size(uk)

    uk1 = selectdim(uk, 1, 1)
    ukN = selectdim(uk, 1, N[1])

    @. uk1 = sqrt(0.5) * uk1
    @. ukN = sqrt(0.5) * ukN

    @inbounds for (i, k) in enumerate(kiter)
        kb = _kbin(k)
        S.sk[kb] = S.sk[kb] + abs2(uk[i])
    end
    norm = 1.e0 / (2 * (N[1] - 1) * prod(N[2:end]))^2
    @. S.sk = norm * S.sk

    if preserve
        @. uk1 = uk1 / sqrt(0.5)
        @. ukN = ukN / sqrt(0.5)
    end
    nothing
end

xfftfreq(nx, ::Val{true})  = rfftfreq(2 * (nx - 1), 2 * (nx - 1))
xfftfreq(nx, ::Val{false}) = fftfreq(nx, nx)

function setup(uk::T; isukreal = true, L = 2π, kwargs...) where {T<:Array{ComplexF64}}
    N = size(uk)
    kiter, kmax = _kgrid(N..., ((nx) -> xfftfreq(nx, Val(isukreal))))
    S = Spectrum(kmax, length(N); L = L)
    return kiter, kmax, S
end

function spectrum(uk::Array{ComplexF64}; preserve = true, kwargs...)
    kiter, kmax, S = setup(uk; kwargs...)
    spectrum!(S, kiter, uk; preserve)
    return S
end

function spectrum(uk::Tuple{Vararg{Array{ComplexF64}}}; preserve = true, kwargs...)
    kiter, kmax, S = setup(uk[1]; kwargs...)
    for uki in uk
        spectrum!(S, kiter, uki; preserve)
    end
    return S
end

function correlation(S::Spectrum; N = 128)
    C = Correlator(N, S.D, typeof(S.Δk); L = S.L/2)
    if D == 2
        @inbounds for (i, r) in enumerate(C.r)
            C.cr[i] = sum((@. S.sk * besselj0(r * S.k)))
        end
        @. C.cr = C.cr / C.cr[1]
    else
        throw(DomainError())
    end
    return C
end
