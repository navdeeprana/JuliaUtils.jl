using IJulia
using PyCall
using PyPlot: plt, Figure, matplotlib

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

function IJulia.metadata(x::Figure)
    w, h = (x.get_figwidth(), x.get_figheight()) .* x.get_dpi()
    return Dict("image/png" => Dict("width" => w / 2, "height" => h / 2))
end

function figax(; nx = 1, ny = 1, h = 3, a = 1.6, kwargs...)
    figsize = (a * h * nx, h * ny)
    fig, ax = plt.subplots(ny, nx; figsize = figsize, kwargs...)
    try
        ax = permutedims(ax, (2, 1))
    catch e
        nothing
    end
    nx * ny > 1 ? ax = hcat(ax...) : nothing
    return fig, ax
end

function noticks!(ax; x = false, y = false)
    ax.xaxis.set_visible(x)
    ax.yaxis.set_visible(y)
    return nothing
end

function fontsize!(ax, fontsize)
    for textobj in ax.findobj(match = (x -> x.__module__ == "matplotlib.text"))
        textobj.set_fontsize(fontsize)
    end
    return nothing
end

function formatter!(axis, low, high)
    f = matplotlib.ticker.ScalarFormatter()
    f.set_powerlimits((low, high))
    axis.set_major_formatter(f)
    return nothing
end
