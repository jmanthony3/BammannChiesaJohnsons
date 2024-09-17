using BammannChiesaJohnsons
using Documenter

DocMeta.setdocmeta!(BammannChiesaJohnsons, :DocTestSetup, :(using BammannChiesaJohnsons); recursive=true)

makedocs(;
    modules=[BammannChiesaJohnsons],
    authors="Joby M. Anthony III, Julian Tse Lop Kun",
    sitename="BammannChiesaJohnsons.jl",
    format=Documenter.HTML(;
        canonical="https://jmanthony3.github.io/BammannChiesaJohnsons.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/BammannChiesaJohnsons.jl",
    devbranch="master",
)
