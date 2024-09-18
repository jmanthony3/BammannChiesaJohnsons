using Pkg; Pkg.precompile()
using Documenter
using BammannChiesaJohnsons

DocMeta.setdocmeta!(BammannChiesaJohnsons, :DocTestSetup, :(using BammannChiesaJohnsons); recursive=true)

makedocs(;
    modules=[BammannChiesaJohnsons],
    authors="Joby M. Anthony III, Julian Tse Lop Kun",
    repo="https://github.com/jmanthony3/BammannChiesaJohnsons.jl/blob/{commit}{path}#{line}",
    sitename="BammannChiesaJohnsons.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/BammannChiesaJohnsons.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/BammannChiesaJohnsons.jl",
    devbranch="main",
)
