# File carries out an entropy analysis, i.e.
# the classical entropy in-equality and Dafermos'
# principle of maximal entropy dissipation is tested
# for a shock calculation

# load solvers
include("sspdg.jl")
include("godsolve.jl")
using CairoMakie

# end time of the calculation
tmax = 1.0
# number of instants in time where the violation of the entropy
# equality is plotted as a heatmap
tsteps = 100
# Number of cells used for the test
K = 100
# number of cells used in the reference Godunov solver
Kgod = K^2
# cfl number for the DG scheme
cfl = 0.02
# amplification factor of the error indicator
fac = 1.0
# calculate DG solutions that should be analysed
grid, sol = NTBurgDGSSP(K, tmax, u1f,fac, cfl)
dgrid, dsol = NTBurgDiscDG(K, tmax, u1f, cfl)

# Matrices to save the entropy
Emat = zeros(K, tsteps)
# and the entropy dissipation
dEmat = zeros(K, tsteps)
# vector to save the total entropy
Etot = zeros(tsteps)
# and the derivative of the total entropy
dEtot = zeros(tsteps)

# set the fontsize of the ploting library for readability.
fontsize_theme = Theme(fontsize = 28, colormap=:grayC)
set_theme!(fontsize_theme)

# loop for the calculation of the total entropies 
# and the violation of the entropy equality
for i=1:tsteps
	dE, dF, E = BurgersdE(grid, sol, (i-1)/tsteps*tmax, fac, cfl)
	Emat[:, i] = dE-dF
	Etot[i] = sum(E)
	ddE, ddF, dE = BurgersdE(grid, dsol, (i-1)/tsteps*tmax, fac, cfl)
        dEmat[:, i] = ddE-ddF
        dEtot[i] = sum(dE)
end

# array of time instants, as a reference for plots
tar = collect(0.0:tmax/tsteps:tmax)[1:tsteps]

# heatmap of the negative logarithm of the violation of the entropy equality
# the small value 1.0E-50 is added for regularization
fig = Makie.heatmap(1:K, tar, log10.(max.(0.0, -Emat) .+ 1.0E-50), axis = (xlabel = "Cells", ylabel="Time"))
Colorbar(fig.figure[1, 2], fig.plot)
CairoMakie.activate!()
save("entroana-"*string(order)*".pdf", fig)

# heatmap of the positive logarithm of the violation of the entropy equlity
fig = Makie.heatmap(1:K, tar, log10.(max.(0.0, Emat) .+ 1.0E-50), axis = (xlabel = "Cells", ylabel="Time"))
Colorbar(fig.figure[1, 2], fig.plot)
CairoMakie.activate!()
save("entroana+"*string(order)*".pdf", fig)

# Array for the total entropy of a fine-grained Godunov reference solution
Etotg = zeros(tsteps)
# calculation of the reference solution
sol = GS(Kgod, tmax, u1f, 2.0, 0.5)
# filling of the array with total entropies
for i=1:tsteps
       Etotg[i] =  sum(U.(sol((i-1)/tsteps*tmax)))/(Kgod)*2
end

# figure with the total entropies over time
fig = lines(tar, Etot, axis = (xlabel = "time", ylabel = "total entropy"), label = "total entropy, DDG", linestyle=:solid, color=:black)
lines!(tar, dEtot, label = "total entropy, DRKDG", linestyle = :dash, color=:black)
lines!(tar, Etotg, label = "Godunov solver total entropy", linestyle=:dot, color=:black)
axislegend(position=:lb)
save("totentroana"*string(order)*".pdf", fig)

GLMakie.activate!()

