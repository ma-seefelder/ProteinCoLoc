using ProteinCoLoc
using Documenter

DocMeta.setdocmeta!(ProteinCoLoc, :DocTestSetup, :(using ProteinCoLoc); recursive=true)

makedocs(;
    modules=[ProteinCoLoc],
    authors="ma-seefelder <manuel.seefelder@uni-ulm.de> and contributors",
    sitename="ProteinCoLoc.jl",
    format=Documenter.HTML(;
        canonical="https://ma-seefelder.github.io/ProteinCoLoc.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ma-seefelder/ProteinCoLoc.jl",
    devbranch="master",
)
