# File runs all tests for the semi-discrete and discrete scheme

# convergence analysis for the semi-discrete scheme
include("convana.jl")
# convergence analysis for the discrete scheme
include("dconvana.jl")
# Pictures of some shock capturing calucaltions to show
# the possibillities to calculate shocks
include("mkntpics.jl")
# the same for the discrete scheme
include("dmkntpics.jl")
# Analysis of the maximal timestep before the calculations of a 
# shock capturing calculation blow up
include("timestepana.jl")
# the same for the discrete scheme
include("dtimestepana.jl")
# test of the semidiscrete entropy inequality
include("entrotest.jl")

