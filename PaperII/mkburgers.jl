# This script carries out the tests for Burgers' equation

# inclusion of the DG solver
include("burgersdg.jl")

# inclusion of a Godunov based reference solver
include("godsolve.jl")

# number of cells for the reference solution
Kref = 30000

rgrid = collect(0.0:2.0/Kref:2.0)[1:Kref] .+ 1/(2*Kref)

# plotting preferences
CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)


# this function carries out a test for Burgers equation
# with the testcase tc, until tend and ktest cells.
# The result is saved under <name>
# the legend can be positioned with pos
function mktestGod(tc, tend, ktest, name; pos=:rt)
        cfl = 0.1/(order^2 + 1)
        grid, sol = NTDGD(BurgersDGDpe!, ktest, tend, tc, cfl=cfl)
        solGod = GS(Kref, tend, tc, 2.0, 0.5)

        fig = PlotSol(grid, sol, tend, 1, yl=L"u(x, t)")
        lines!(rgrid, solGod(tend)[:, 1], label = "reference", linestyle=:dash)
        axislegend(position=pos)
        save("figures/"*name*string(ktest)*"Burgers"*string(order)*".pdf", fig)

end

# Sine initial condition
mktestGod(u1f, 0.5, 20, "sin05")
mktestGod(u1f, 1.0, 20, "sin10", pos = :lt)
mktestGod(u1f, 0.5, 50, "sin05")
mktestGod(u1f, 1.0, 50, "sin10", pos = :lt)



# Sonic point IC
mktestGod(u2f, 0.3, 20, "sp03", pos =:lt)
mktestGod(u2f, 0.3, 50, "sp03", pos =:lt)

