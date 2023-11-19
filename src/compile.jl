using PackageCompiler
create_app(
    "C:/Users/manue/Documents/GitHub/ProteinCoLoc", 
    "C:/Users/manue/Desktop/ProteinCoLoc", 
    script="C:/Users/manue/Documents/GitHub/ProteinCoLoc/src/script.jl";
    #executables = ["ProteinCoLoc" => "julia_main"],
    force = true
    )