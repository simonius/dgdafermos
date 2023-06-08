# File produces a set of example solutions for the shock capturing 
# capability of the method
# include implementation of the semi-discrete DG solver
include("eulerdg.jl")

# include eno2 solver for reference purposes
include("lfsolver.jl")

# Number of cells in the reference solution
Kref = 30000

# reference grid
rgrid = collect(0.0:10.0/Kref:10.0)[1:Kref] .+ 1/(2*Kref)

CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)


function mktest(tc, tend, ktest, name, cfl = 0.1/(order^2 + 1))
	grid, sol = NTDGD(EulerDGDpe!, ktest, tend, tc, cfl=cfl)
	solLF = EulerLFlu(tc, tend, Kref, 5.0)

	fig = PlotSol(grid, sol, tend, 1)
	lines!(rgrid, solLF[:, 1], label = "reference", linestyle=:dash)
	axislegend()
	save("figures/"*name*"1k"*string(ktest)*"p"*string(order)*".pdf", fig)

	fig=PlotSol(grid, sol, tend, 2)
	lines!(rgrid, solLF[:, 2], label = "reference", linestyle=:dash)
	axislegend()
	save("figures/"*name*"2k"*string(ktest)*"p"*string(order)*".pdf", fig)

	fig=PlotSol(grid, sol, tend, 3)
	lines!(rgrid, solLF[:, 3], label = "reference", linestyle=:dash)
	axislegend()
	save("figures/"*name*"3k"*string(ktest)*"p"*string(order)*".pdf", fig)
end



######  First test: shock tube with 33 and 100 cells.
mktest(shuosh6a, 2.0, 13, "shuosh6a")
mktest(shuosh6a, 2.0, 25, "shuosh6a")
mktest(shuosh6a, 2.0, 100, "shuosh6a")

######  Second test: shock tube with 33 and 100 cells.
mktest(shuosh6b, 1.4, 13, "shuosh6b")
mktest(shuosh6b, 1.4, 25, "shuosh6b")
mktest(shuosh6b, 1.4, 100, "shuosh6b")

######  Third test: Shu-Osher with 50, 100 and 200 cells.

mktest(shuosh6, 1.8, 50, "shuosh6")
mktest(shuosh6, 1.8, 100, "shuosh6")

if order < 5
	mktest(shuosh6, 1.8, 200, "shuosh6")
else
	mktest(shuosh6, 1.8, 25, "shuosh6")
end
