using PackageCompiler
create_app(
    "C:/Users/manue/Documents/GitHub/ProteinCoLoc", 
    "C:/Users/manue/Desktop/ProteinCoLoc", 
    script="script.jl";
    executables = ["ProteinCoLoc" => "julia_main"]
    )