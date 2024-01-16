using LaTeXStrings


# File produces a set of example solutions for the shock capturing 
# capability of the method for the Euler equations

# include implementation of the semi-discrete DG solver
include("eulerdg.jl")

# include eno2 solver for reference purposes
include("lfsolver.jl")

# exact Rieman solver for reference purposes
include("rsolver.jl")

# Number of cells in the reference solution
Kref = 30000

# reference grid
rgrid = collect(0.0:10.0/Kref:10.0)[1:Kref] .+ 1/(2*Kref)

# plotting preferences
CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

# tests the DG solver for the Euler equations with a reference solution based on a MLF solver
function mktestLF(tc, tend, ktest, name; cfl = 0.1/(order^2 + 1), lpr=:rt, lpm=:rt, lpe=:rt)
	grid, sol = NTDGD(EulerDGDpe!, ktest, tend, tc, cfl=cfl)
	solLF = EulerLFlu(tc, tend, Kref)

	fig = PlotSol(grid, sol, tend, 1, yl=L"\rho(x, t)")
	lines!(rgrid, solLF[:, 1], label = "reference", linestyle=:dash)
	axislegend(position=lpr)
	save("figures/"*name*"1k"*string(ktest)*"p"*string(order)*".pdf", fig)

	fig=PlotSol(grid, sol, tend, 2, yl=L"\rho(x, t) v(x, t)")
	lines!(rgrid, solLF[:, 2], label = "reference", linestyle=:dash)
	axislegend(position=lpm)
	save("figures/"*name*"2k"*string(ktest)*"p"*string(order)*".pdf", fig)

	fig=PlotSol(grid, sol, tend, 3, yl=L"E(x, t)")
	lines!(rgrid, solLF[:, 3], label = "reference", linestyle=:dash)
	axislegend(position=lpe)
	save("figures/"*name*"3k"*string(ktest)*"p"*string(order)*".pdf", fig)
end

# tests the DG solver for the Euler equations with a reference solution provided by an exact Riemann
# solver (only usable for Riemann initial data)
function mktestRS(tc, tend, ktest, name, x0; cfl = 0.1/(order^2 + 1), lpr=:rt, lpm=:rt, lpe=:rt)
        grid, sol = NTDGD(EulerDGDpe!, ktest, tend, tc, cfl=cfl)
        ul = tc(0.0)
        ur = tc(10.0)
        solRS1(x) = uRiemann(ul, ur, (x - x0)/tend)[1]
        solRS2(x) = uRiemann(ul, ur, (x - x0)/tend)[2]
        solRS3(x) = uRiemann(ul, ur, (x - x0)/tend)[3]

        rg = collect(0.0:0.01:10.0)


        fig = PlotSol(grid, sol, tend, 1,  yl=L"\rho(x, t)")
        lines!(rg, solRS1.(rg), label = "reference", linestyle=:dash)
        axislegend(position=lpr)
        save("figures/"*name*"1k"*string(ktest)*"p"*string(order)*".pdf", fig)

        fig=PlotSol(grid, sol, tend, 2, yl=L"\rho(x, t) v(x, t)")
        lines!(rg, solRS2.(rg), label = "reference", linestyle=:dash)
        axislegend(position=lpm)
        save("figures/"*name*"2k"*string(ktest)*"p"*string(order)*".pdf", fig)

        fig=PlotSol(grid, sol, tend, 3, yl=L"E(x, t)")
        lines!(rg, solRS3.(rg), label = "reference", linestyle=:dash)
        axislegend(position=lpe)
        save("figures/"*name*"3k"*string(ktest)*"p"*string(order)*".pdf", fig)
end



######  First test: Sod shock tube with 13, 25 and 100 cells.
mktestRS(shuosh6a, 2.0, 13, "shuosh6a", 5.0, lpm=:lt)
mktestRS(shuosh6a, 2.0, 25, "shuosh6a", 5.0, lpm=:lt)
mktestRS(shuosh6a, 2.0, 100, "shuosh6a", 5.0, lpm=:lt)

######  Second test: Toro shock tube 1
mktestRS(toro1, 2.0, 13, "toro1", 3.0)
mktestRS(toro1, 2.0, 25, "toro1", 3.0)
mktestRS(toro1, 2.0, 100, "toro1", 3.0)

######  Third test: Lax shock tube with 13, 25 and 100 cells.
mktestRS(shuosh6b, 1.4, 13, "shuosh6b", 5.0, lpr=:lt, lpm=:lt, lpe=:lb)
mktestRS(shuosh6b, 1.4, 25, "shuosh6b", 5.0, lpr=:lt, lpm=:lt, lpe=:lb)
mktestRS(shuosh6b, 1.4, 100, "shuosh6b", 5.0, lpr=:lt, lpm=:lt, lpe=:lb)

###### 4-th test: Toro shock tube 3
mktestRS(toro3, 0.12, 13, "toro3", 5.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:lt, lpe=:lb)
mktestRS(toro3, 0.12, 25, "toro3", 5.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:lt, lpe=:lb)
mktestRS(toro3, 0.12, 100, "toro3", 5.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:lt, lpe=:lb)

###### 5-th test: Toro shock tube 4
mktestRS(toro4, 0.35, 13, "toro4", 4.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:lt, lpe=:lt)
mktestRS(toro4, 0.35, 25, "toro4", 4.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:lt, lpe=:lt)
mktestRS(toro4, 0.35, 100, "toro4", 4.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:lt, lpe=:lt)

###### 6-th test: Toro shock tube 5
mktestRS(toro5, 0.12, 13, "toro5", 8.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:cb)
mktestRS(toro5, 0.12, 25, "toro5", 8.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:cb)
mktestRS(toro5, 0.12, 100, "toro5", 8.0, cfl=0.005/(order^2+1), lpr=:lt, lpm=:cb)


######  Last test: Shu-Osher with 50, 100 and 200 cells.

mktestLF(shuosh6, 1.8, 50, "shuosh6", lpr=:lb)
mktestLF(shuosh6, 1.8, 100, "shuosh6", lpr=:lb)


mktestLF(shuosh6, 1.8, 200, "shuosh6", lpr=:lb)

mktestLF(shuosh6, 1.8, 25, "shuosh6", lpr=:lb)

