# inclusion of this script carries out the tests for the Buckley-Leverett equation

# inclusion of the DG solver for the BL equation
include("bldg.jl")

# The number of cells for the reference solution calculation
Kref = 30000

# setting the plotting preferences
CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

# Solver used to produce the reference solution plots
# based on a first order HLL solver. 
function BLHLL!(du, u, p, t)
	dx = p[1]
	N = length(du)
	for n=1:length(du)
		du[n] = (hBLHLL(u[mod1(n-1, N)], u[n]) - hBLHLL(u[n], u[mod1(n+1, N)]))/dx
	end
end

# driver routine for the reference solution solver.
# Kref is the number of cells,
# tend the end time of the simulation,
# tc the initial condition
# xmax the right end of the domain,
# cfl the ratio between dt and dx
function BLref(Kref, tend, tc, xmax, cfl=0.1)
	rgrid = collect(0.0:xmax/Kref:xmax)[1:Kref] .+ xmax/(2*Kref)
	u0 = zeros(Kref)
	for n=1:Kref
		u0[n] = tc(rgrid[n])[1]
	end
	dx = rgrid[2] - rgrid[1]
	dt = cfl*dx
	p = [dx]
	prob = ODEProblem(BLHLL!, u0, (0.0, tend), p)
	sol = solve(prob, Euler(), dt =dt, adaptive=false, saveat=[0.0, tend])
	return rgrid, sol
end

# this routine tests the DG solver for the Buckley-Leverett equation against
# a HLL based solver. 
# tc is the testcase
# tend the time at which the results are compared,
# ktest is the number of cells used for the DG solver
# the resulting images are saved with the filename <name>
function mktestBL(tc, tend, ktest, name)
        cfl = 0.1/(order^2 + 1)
        grid, sol = NTDGD(BLDGDpe!, ktest, tend, tc, cfl=cfl)
        rgrid, solHLL = BLref(Kref, tend, tc, 2.0)

        fig = PlotSol(grid, sol, tend, 1, yl=L"u(x, t)")
        lines!(rgrid, solHLL(tend), label = "reference", linestyle=:dash)
        axislegend(position=:lt)
        save("figures/"*name*string(ktest)*"BL"*string(order)*".pdf", fig)
end

# calls to the function above to test the solver
mktestBL(u1f, 0.6, 13, "BLRiemann")
mktestBL(u1f, 0.6, 25, "BLRiemann")
mktestBL(u1f, 0.6, 50, "BLRiemann")

