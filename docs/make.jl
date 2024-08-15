# set up environment
import Pkg; Pkg.activate("./docs")
using Documenter, ProteinCoLoc
###############################################################################
DocMeta.setdocmeta!(ProteinCoLoc, :DocTestSetup, :(using ProteinCoLoc); recursive=true)

makedocs(;
    modules=[ProteinCoLoc],
    authors="ma-seefelder <manuel.seefelder@uni-ulm.de>",
    sitename="ProteinCoLoc.jl",
    format=Documenter.HTML(;
        canonical="https://ma-seefelder.github.io/ProteinCoLoc.jl",
        edit_link="master",
        assets=String["assets/custom_css.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Basic usage" => "basic_usage.md",
        "API" => [
            "Main function" => "main.md",
            "Load Images" => "loading_image.md",
            "Plot" => "plot.md",
            "Infer colocalization" => "colocalisation.md"
        ], 
    ],
    checkdocs = :none,
    pagesonly = true,
    remotes = nothing
)

#deploydocs(;repo="github.com/ma-seefelder/ProteinCoLoc")
