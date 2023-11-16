function julia_main()::CInt
    include("ProteinCoLoc.jl")
    using .ProteinCoLoc
    ProteinCoLoc.gui()
    return 0
end