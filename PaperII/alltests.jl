# This script carries out all implemented tests for the Euler equations
include("eulerdg.jl")

# convergence analysis for the semi-discrete scheme
include("convana.jl")
# Pictures of some shock capturing calucaltions to show
# the possibillities to calculate shocks
include("mkntpics.jl")
# Analysis of the maximal timestep before the calculations of a 
# shock capturing calculation blows up
include("timestepana.jl")
# test of the semidiscrete entropy inequality
include("entrotest.jl")

