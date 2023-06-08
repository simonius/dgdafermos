# File produces a set of example solutions for the shock capturing 
# capability of the method
# include implementation of the semi-discrete DG solver
include("sspdg.jl")

# include Godunov solver for reference purposes
include("godsolve.jl")

# Number of cells in the reference solution
Kref = 1000 
# Number of cells for two test runs
Ktest1 = 20
Ktest2 = 50

# Factor of misprediction for test purposes
# set to 1.0 for the solver as reported in the publication
mfac = 1.0

# reference grid
rgrid = 0.0:2/Kref:2.0-2/Kref

CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)


######  First tests: sine at 20
grid, sol = NTBurgDGSSP(Ktest1, 1.0, u1f, mfac)
rsol = GS(Kref, 1.0, u1f, 2.0, 0.5)


fig = PlotSol(grid, sol, 0.25, ls=:dot)
lines!(rgrid, rsol(0.25), label = "reference")
axislegend()
save("sin1at025p"*string(order)*".pdf", fig)

fig=PlotSol(grid, sol, 0.5, ls=:dot)
lines!(rgrid, rsol(0.5), label = "reference")
axislegend()
save("sin1at050p"*string(order)*".pdf", fig)

fig=PlotSol(grid, sol, 1.0, ls=:dot)
lines!(rgrid, rsol(1.0), label = "reference")
axislegend()
save("sin1at100p"*string(order)*".pdf", fig)

######  Second tests: sine at ktest2

grid, sol = NTBurgDGSSP(Ktest2, 1.0, u1f, mfac)
fig = PlotSol(grid, sol, 0.25, ls=:dot)
lines!(rgrid, rsol(0.25), label = "reference")
axislegend()
save("sin2at025p"*string(order)*".pdf", fig)

fig = PlotSol(grid, sol, 0.5, ls=:dot)
lines!(rgrid, rsol(0.5), label = "reference")
axislegend()
save("sin2at050p"*string(order)*".pdf", fig)

fig = PlotSol(grid, sol, 1.0, ls=:dot)
lines!(rgrid, rsol(1.0), label = "reference")
axislegend()
save("sin2at100p"*string(order)*".pdf", fig)


######  Sonic point glitch tests
grid, sol = NTBurgDGSSP(Ktest1, 1.0, u2f, mfac)
rsol = GS(Kref, 1.0, u2f, 2.0, 0.25)


fig = PlotSol(grid, sol, 0.3, ls=:dot)
lines!(rgrid, rsol(0.3), label = "reference")
axislegend(position=:rb)
save("rare1at030p"*string(order)*".pdf", fig)


grid, sol = NTBurgDGSSP(Ktest2, 1.0, u2f, mfac)
fig = PlotSol(grid, sol, 0.3, ls=:dot)
lines!(rgrid, rsol(0.3), label = "reference")
axislegend(position=:rb)
save("rare2at030p"*string(order)*".pdf", fig)




