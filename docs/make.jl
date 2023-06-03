using JuliaUtils
using Documenter

DocMeta.setdocmeta!(JuliaUtils, :DocTestSetup, :(using JuliaUtils); recursive=true)

makedocs(;
    modules=[JuliaUtils],
    authors="Navdeep Rana",
    repo="https://github.com/navdeeprana/JuliaUtils.jl/blob/{commit}{path}#{line}",
    sitename="JuliaUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://navdeeprana.github.io/JuliaUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
    "Home" => "index.md",
	"Core utility functions" => [ 
        "I/O"=> "core/io.md",
        "Space"=> "core/space.md",
        "Periodic"=> "core/periodic.md",
        "Physics"=> "core/physics.md",
        "Defects"=> "core/defects.md",
        "Spectrum"=> "core/spectrum.md",
        "Interpolations"=> "core/interpolations.md"
        ],
	"Ploting Utilities" => [
        "Makie"    => "ploting/makie.md",
        "PyPlot"   => "ploting/pyplot.md"
        ],
	"Macros" =>  "macros/macros.md",
    "Library" => [
        "Public" => "library/public.md", 
        "Internals" => "library/internals.md",
        "Glossary" => "library/glossary.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/navdeeprana/JuliaUtils.jl",
    devbranch="main",
)
