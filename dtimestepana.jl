# This file tests the discrete scheme for the maximal allowable time-step
# Simulations with time-steps and cells numbers on a regular grid are
# started and the maximal simulation time is reported
include("sspdg.jl")
using Makie
using CairoMakie
using LaTeXStrings

# maximal simulation time after which the simulation is stopped,
# even when no blow-up occured
tmax = 1.0

# minimum factor between time- and space-step, multiplied by $p^2 + 1$
CFLstart = 0.2
# increase
CFLstep = 0.2
# maximum factor to test
CFLend = 4.0

CFLarray = (CFLstart:CFLstep:CFLend)/(order^2 + 1)

# Array of tested amounts of cells.
Karray = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]

n = length(CFLarray)
m = length(Karray)

# array holding the blow-up time of the simulation
tendar = zeros(n, m)

# Loop tests the grid of cell numbers and time-steps
# saves the final simulation time in the tendar array
for i=1:n
        for j=1:m
                print(Karray[j], " ", CFLarray[i])
                grid, sol = NTBurgDiscDG(Karray[j], tmax, u1f, CFLarray[i])
                tendar[i, j] = sol.t[end]   
        end
end

# set font size for readability
fontsize_theme = Theme(fontsize = 28, colormap=:grayC)
set_theme!(fontsize_theme)
# heatmap of the end times
fig = Makie.heatmap(CFLarray.*(order^2+1), Karray, tendar, axis = (xlabel = L"\frac{\Delta t}{\Delta x}(p^2 + 1)", ylabel="K"))
Colorbar(fig.figure[1, 2], fig.plot)
CairoMakie.activate!()
save("dcflana"*string(order)*".pdf", fig) 
