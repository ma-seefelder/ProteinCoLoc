#using BenchmarkTools
using .ProteinCoLoc

@time ProteinCoLoc.start_analysis(
    "C:/Users/manue/Desktop/benchmark/positive",
    "C:/Users/manue/Desktop/benchmark/negative",
    "C:/Users/manue/Desktop/benchmark/result/", # path to the output folder
    64, # number of patches
    200, # number of patches for local correlation
    3, # number of channels
    true, # channel selection
    [2,3], # channel selection two
    false, # patched correlation plot
    false, # local correlation plot
    false, # bayes factor plot
    false, # bayes range plot
    false, # posterior plot
    false, # mask plot
)
