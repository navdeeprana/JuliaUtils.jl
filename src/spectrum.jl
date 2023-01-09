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

# Returns the maximum possible wavenumber for the range kiter.
# (r)fftfreq omits the largest possible wavenumber so we add one.
@inline function _kmax(kiter)
    kmax = maximum(kiter) .+ 1
    return round(Int, hypot(kmax...)) + 1
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

# Computes the power spectrum S(k) given the fourier transform uk = ⨏(u) of the data.
# TODO: Write documentation.
# The internals differ slightly for a real u and a complex u.
function spectrum!(
    S::Spectrum,
    kiter::Iterators.ProductIterator,
    uk::Array{ComplexF64},
    norm;
    isreal = true,
    preserve = true
)
    N   = size(uk)
    uk0 = selectdim(uk, 1, 1)
    ukN = selectdim(uk, 1, N[1])
    if isreal
        @. uk0 = sqrt(0.5) * uk0
        @. ukN = sqrt(0.5) * ukN
    end

    @inbounds for (i, k) in enumerate(kiter)
        kb = _kbin(k)
        S.sk[kb] = S.sk[kb] + abs2(uk[i])
    end
    @. S.sk = norm * S.sk

    if isreal & preserve
        @. ukN = ukN / sqrt(0.5)
        @. uk0 = uk0 / sqrt(0.5)
    end
    nothing
end

xfftfreq(nx, ::Val{true})  = rfftfreq(2 * (nx - 1), 2 * (nx - 1))
xfftfreq(nx, ::Val{false}) = fftfreq(nx, nx)

function setup(uk::Array{ComplexF64}; isreal = true, L = 2π)
    N = size(uk)
    kiter, kmax = _kgrid(N..., ((nx) -> xfftfreq(nx, Val(isreal))))
    if isreal
        norm = 1.0 / (2 * (N[1] - 1) * prod(N[2:end]))^2
    else
        norm = 1.0 / (2.0 * prod(N)^2)
    end
    S = Spectrum(kmax, length(N); L = L)
    return kiter, kmax, norm, S
end

function spectrum(uk::Array{ComplexF64}; preserve = true, isreal = true, kwargs...)
    kiter, kmax, norm, S = setup(uk; isreal, kwargs...)
    spectrum!(S, kiter, uk, norm; isreal, preserve)
    return S
end

function spectrum(uk::Tuple{Vararg{Array{ComplexF64}}}; preserve = true, isreal = true, kwargs...)
    kiter, kmax, norm, S = setup(uk[1]; isreal, kwargs...)
    for uki in uk
        spectrum!(S, kiter, uki, norm; isreal, preserve)
    end
    return S
end

function correlation(S::Spectrum; N = 128)
    C = Correlator(N, S.D, typeof(S.Δk); L = S.L / 2)
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
