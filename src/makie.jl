using Makie: Figure, Axis, AxisAspect

function figax(; nx = 1, ny = 1, h = 3, a = 1.6, s = 100, kwargs...)
    res = (round(Int, a * s * h * nx), round(Int, s * h * ny))
    fig = Figure(resolution = res)
    ax = [Axis(fig[j, i]; aspect = AxisAspect(a), kwargs...) for i in 1:nx, j in 1:ny]
    if nx * ny == 1
        return fig, ax[1]
    else
        return fig, ax
    end
end
