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
    ],
)

deploydocs(;
    repo="github.com/navdeeprana/JuliaUtils.jl",
    devbranch="main",
)
