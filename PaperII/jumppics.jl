# This script is used to produce the Figure 4 in the revised manuscript
# Figure 4 showcases that the L^2 projection of a discontinouos solution on a space
# of piecewise polynomials results in wrong jump hights in the polynomial 
# representation. This effect gets worse for growing orders.

# inclusion of the basis functions
include("bases.jl")
using Makie
using CairoMakie
using LaTeXStrings

# preferences for ploting
CairoMakie.activate!()
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)


# the projection routine projects the function uf onto the
# first Nmax Legendre polynomials
function L2Proj(uf, Nmax)
        c = zeros(Nmax)
        for k=1:Nmax
                integ(x) = uf(x)*Pleg(k-1, x)
                int, err = quadgk(integ, -1.0, 1.0, maxevals=100)
                c[k] = (2*k-1)*int*0.5
        end
        return c
end


# Plots the Ansatzfunctions of 4 elements
# when the function uf is projected onto them.
function jumppics(uf)
	K = 4
	xl = -1.0
	xr = 1.0
	dx = (xr - xl)/K
	grid = xl .+ dx*collect(1:K) .- dx/2
	order = 7
	minord = 4
	c = zeros(K, order+1)
        for n=1:K
                c[n, :] = L2Proj(x->uf(grid[n] + dx/2*x), order+1)
	end

	rgrid = collect(-1.0:0.001:1.0)
	axis = (xlabel = L"x", ylabel = L"y")
	styles = reverse([:dot, :dash, :dashdot, :solid])
	p = scatter(rgrid, uf.(rgrid), color=:black, label = L"u(x)", axis = axis)
	for ord=minord:order
		for k=1:K
			locgrid = grid[k] .+ dx/2*rgrid
			pf(x) = CellEval(c[k, 1:ord], x)
			if k == 1
				o = ord
				lines!(locgrid, pf.(rgrid), color=:black, label = L"\mathrm{P}_{%$o} u(x)", linestyle=styles[ord-minord+1] )
			else
				lines!(locgrid, pf.(rgrid), color=:black, linestyle=styles[ord-minord+1])
			end
		end
	end

	axislegend(position=:lt)
	return p
end

# the used example function. Please note that the jump is located 
# 1/100th parts right of the cell boundary
uf(x) = x < 0.01 ? 0.0 : 1.0

fig = jumppics(uf)
save("figures/jumpfig.pdf", fig)
