# This file tests the semi-discrete scheme for the maximal allowable time-step
# Simulations with time-steps and cells numbers on a regular grid are
# started and the maximal simulation time is reported

include("eulerdg.jl")
using Makie
using CairoMakie
using LaTeXStrings

# maximal simulation time after which the simulation is stopped,
# even when no blow-up occured
tmax = 1.0

# minimum factor between time- and space-step, multiplied by $p^2 + 1$
CFLstart = 0.1
# increase
CFLstep = 0.1
# maximum factor to test
CFLend = 2.0

CFLarray = (CFLstart:CFLstep:CFLend)/(order^2 + 1)

# Array of tested amounts of cells
Karray = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]


n = length(CFLarray)
m = length(Karray)

# array holding blow-up times
tendar = zeros(n, m)

# Loop tests the grid of cell numbers and time-steps
# saves the final simulation time in the tendar array
for i=1:n
	println("cfl = ", CFLarray[i])
	print("K = ")

        Threads.@threads for j=1:m
		print(Karray[j],  " ")
                grid, sol = NTDGD(EulerDGDpe!, Karray[j], tmax, shuosh6a, cfl = CFLarray[i])
                tendar[i, j] = sol.t[end]   
        end
end


fontsize_theme = Theme(fontsize = 28, colormap=:grayC)
set_theme!(fontsize_theme)

# heatmap of end times
fig = Makie.heatmap(CFLarray.*(order^2+1), Karray, tendar, axis = (xlabel = L"\frac{\Delta t}{\Delta x}(p^2 + 1)", ylabel="K"))
Colorbar(fig.figure[1, 2], fig.plot)
CairoMakie.activate!()
save("figures/cflana"*string(order)*".pdf", fig) 
