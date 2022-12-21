# File implements the printing of same example
# shock calculations using the discrete scheme

# load implementation of the DG and reference Godunov solver
include("sspdg.jl")
include("godsolve.jl")

# Number of cells used in the reference solution
Kref = 1000 

# Number of cells used in the two test runs
Ktest1 = 20
Ktest2 = 50

rgrid = 0.0:2/Kref:2.0-2/Kref

CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)


######  First tests: sine at 20
grid, sol = NTBurgDiscDG(Ktest1, 1.0, u1f)
rsol = GS(Kref, 1.0, u1f, 2.0, 0.5)


fig = PlotSol(grid, sol, 0.25, ls=:dot)

lines!(rgrid, rsol(0.25), label = "reference")
axislegend()
save("dsin1at025p"*string(order)*".pdf", fig)

fig=PlotSol(grid, sol, 0.5, ls=:dot)
lines!(rgrid, rsol(0.5), label = "reference")
axislegend()
save("dsin1at050p"*string(order)*".pdf", fig)

fig=PlotSol(grid, sol, 1.0, ls=:dot)
lines!(rgrid, rsol(1.0), label = "reference")
axislegend()
save("dsin1at100p"*string(order)*".pdf", fig)

######  Second tests: sine at ktest2

grid, sol = NTBurgDiscDG(Ktest2, 1.0, u1f)
fig = PlotSol(grid, sol, 0.25, ls=:dot)
lines!(rgrid, rsol(0.25), label = "reference")
axislegend()
save("dsin2at025p"*string(order)*".pdf", fig)

fig = PlotSol(grid, sol, 0.5, ls=:dot)
lines!(rgrid, rsol(0.5), label = "reference")
axislegend()
save("dsin2at050p"*string(order)*".pdf", fig)

fig = PlotSol(grid, sol, 1.0, ls=:dot)
lines!(rgrid, rsol(1.0), label = "reference")
axislegend()
save("dsin2at100p"*string(order)*".pdf", fig)


######  Sonic point glitch tests
grid, sol = NTBurgDiscDG(Ktest1, 1.0, u2f)
rsol = GS(Kref, 1.0, u2f, 2.0, 0.25)


fig = PlotSol(grid, sol, 0.3, ls=:dot)
lines!(rgrid, rsol(0.3), label = "reference")
axislegend()
save("drare1at030p"*string(order)*".pdf", fig)


grid, sol = NTBurgDiscDG(Ktest2, 1.0, u2f)
fig = PlotSol(grid, sol, 0.3, ls=:dot)
lines!(rgrid, rsol(0.3), label = "reference")
axislegend()
save("drare2at030p"*string(order)*".pdf", fig)




