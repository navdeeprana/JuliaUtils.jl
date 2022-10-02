using FFTW
using SpecialFunctions
using Statistics

abstract type AbstractCorrelator{T} end
abstract type AbstractSpectrum{T} end

struct Correlator{T<:AbstractFloat} <: AbstractCorrelator{T}
    D::Integer
    r::Vector{T}
    cr::Vector{T}
end

function Correlator(N, D, T; L = π)
    r = LinRange(0, L, N + 1)
    cr = zeros(T, N + 1)
    return Correlator{T}(D, r, cr)
end

struct Spectrum{T<:AbstractFloat} <: AbstractSpectrum{T}
    D::Integer
    L::T
    dk::T
    k::Vector{T}
    sk::Vector{T}
end

function Spectrum(N, D; L = 2π)
    dk = 2π / L
    T = typeof(dk)
    k = @. dk * LinRange(0, N, N + 1)
    sk = zeros(T, N + 1)
    return Spectrum{T}(D, L, dk, k, sk)
end

function kgrid(nx::T; hermitian = true) where {T<:Integer}
    if hermitian
        kiter = Iterators.product(rfftfreq(2 * (nx - 1), 2 * (nx - 1)))
    else
        kiter = Iterators.product(fftfreq(nx, nx))
    end
    return kiter, round(Int, hypot(maximum(kiter)...)) + 1
end

function kgrid(nx::T, ny::T; hermitian = true) where {T<:Integer}
    if hermitian
        kiter = Iterators.product(rfftfreq(2 * (nx - 1), 2 * (nx - 1)), fftfreq(ny, ny))
    else
        kiter = Iterators.product(fftfreq(nx, nx), fftfreq(ny, ny))
    end
    return kiter, round(Int, hypot(maximum(kiter)...)) + 1
end

function kgrid(nx::T, ny::T, nz::T; hermitian = true) where {T<:Integer}
    if hermitian
        kiter = Iterators.product(
            rfftfreq(2 * (nx - 1), 2 * (nx - 1)),
            fftfreq(ny, ny),
            fftfreq(nz, nz)
        )
    else
        kiter = Iterators.product(fftfreq(nx, nx), fftfreq(ny, ny), fftfreq(nz, nz))
    end
    return kiter, round(Int, hypot(maximum(kiter)...)) + 1
end

@inline function get_kbin(k::Tuple{Vararg{T}}) where {T<:Real}
    return round(Int, hypot(k...)) + 1
end

function spectrum!(
    S::Spectrum,
    kiter::Iterators.ProductIterator,
    uk::Array{ComplexF64};
    preserve = true
)
    N = size(uk)

    uk1 = selectdim(uk, 1, 1)
    ukN = selectdim(uk, 1, N[1])

    @. uk1 = sqrt(0.5e0) * uk1
    @. ukN = sqrt(0.5e0) * ukN

    @inbounds for (i, k) in enumerate(kiter)
        kbin = get_kbin(k)
        S.sk[kbin] = S.sk[kbin] + abs2(uk[i])
    end
    norm = 1.e0 / (2 * (N[1] - 1) * prod(N[2:end]))^2
    @. S.sk = norm * S.sk

    if preserve
        @. uk1 = uk1 / sqrt(0.5e0)
        @. ukN = ukN / sqrt(0.5e0)
    end
    nothing
end

function spectrum(uk::T; L = 2π, kwargs...) where {T<:Array{ComplexF64}}
    N = size(uk)
    kiter, kmax = kgrid(N...)
    S = Spectrum(kmax, length(N); L = L)
    spectrum!(S, kiter, uk, kwargs...)
    return S
end

function spectrum(u::Tuple{Vararg{T}}; L = 2π, kwargs...) where {T<:Array{ComplexF64}}
    N = size(u[1])
    kiter, kmax = kgrid(N...)
    S = Spectrum(kmax, length(N); L = L)
    for ui in u
        spectrum!(S, kiter, ui, kwargs...)
    end
    return S
end

function correlation(S::Spectrum)
    C = Correlator(128, S.D, typeof(S.dk); L = 0.5e0 * S.L)
    for (i, r) in enumerate(C.r)
        C.cr[i] = sum((@. S.sk * besselj0(r * S.k)))
    end
    @. C.cr = C.cr / C.cr[1]
    return C
end
