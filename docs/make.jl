using ProteinCoLoc
using Documenter

DocMeta.setdocmeta!(ProteinCoLoc, :DocTestSetup, :(using ProteinCoLoc); recursive=true)

makedocs(;
    modules=[ProteinCoLoc],
    authors="ma-seefelder <manuel.seefelder@uni-ulm.de>",
    sitename="ProteinCoLoc.jl",
    format=Documenter.HTsML(;
        canonical="https://ma-seefelder.github.io/ProteinCoLoc.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Load Images" => "loading_image.md"
    ],
    checkdocs = :none
)

#=
deploydocs(;
    repo="github.com/ma-seefelder/ProteinCoLoc",
    devbranch="registration",
)
=#