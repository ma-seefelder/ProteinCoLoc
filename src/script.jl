using .ProteinCoLoc
function julia_main()::CInt
    include("ProteinCoLoc.jl")
    ProteinCoLoc.gui()
end