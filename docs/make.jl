using ProteinCoLoc
using Documenter

DocMeta.setdocmeta!(ProteinCoLoc, :DocTestSetup, :(using ProteinCoLoc); recursive=true)

makedocs(;
    modules=[ProteinCoLoc],
    authors="ma-seefelder <manuel.seefelder@uni-ulm.de>",
    sitename="ProteinCoLoc.jl",
    format=Documenter.HTML(;
        canonical="https://ma-seefelder.github.io/ProteinCoLoc.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Load Images" => "loading_image.md",
        "Plot" => "plot.md"
    ],
    checkdocs = :none
)

#=
deploydocs(;
    repo="github.com/ma-seefelder/ProteinCoLoc",
    devbranch="registration",
)
=#