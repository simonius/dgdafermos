# includes the implementation of the solvers
include("eulerdg.jl")

# File performs a convergence analysis for the semi-discrete scheme
# The initial condition smootheuler, a density variation rotating around 
# the periodic BCs is used
# the reference is calculated using characteristic backtracing

# calculates the Lagrange Polynomial k to x_k at x
function Plag(xar, k, x)
	K = length(xar)
	v = 1.0
	for i=1:k-1
		v = v*(x-xar[i])/(xar[k] - xar[i])
	end 
	for i=k+1:K
		v = v*(x-xar[i])/(xar[k] - xar[i])
	end 
	return v
end

# calcualtes the Value of $u$ at $x$
function SCeval(xar, uT, x)
	v = 0.0
	K = length(xar)
	for k=1:K
		v = v + Plag(xar, k, x)*uT[k]
	end 
	return v
end

# evaluates the solution at x
function SolEval(grid, sol, x)
        solval = zeros(3)
	#find cell
	starts = grid[:, 1]
	n = length(starts)
	lb = 1
	ub = n+1
	im = 1
	while ub - lb > 1
		im = floor(Int, 0.5*(ub + lb))
		if starts[im] < x
			lb = im
		else
			ub = im
		end
	#	println(lb, " ", im, " ", ub)
	end
	# return the interpolated value in cell lb
	for l=1:3
                solval[l] = SCeval(grid[lb, :], sol[lb, :, l], x)
        end
        return solval
end

# the initial condition is a smooth density variation
uinit = smootheuler
# that is advected to the right
refsol(x, t) = uinit(mod(x - 2.0*t, 10.0))


# Points used for Monte-Carlo quadrature of solution error
Nmc = 10000

# End time
tend = 5.0

# Used amount of cells.
if order < 6
	Narr = collect(10:5:100)
else
	Narr = collect(4:2:28)
end

# Used CFL number for the calculation with the lowest number
# of cells, scaled with the number. This is needed as we are scaling 
# down the CFL number together with the cells to allow up to sixt order convergence
cflp = 0.02*Narr[1]

# Arrays to save the errors in the two and one norms.
errarr = zeros(length(Narr))
errarr2 = zeros(length(Narr))
errarrinf = zeros(length(Narr))
errarrUS = zeros(length(Narr))
errarr2US = zeros(length(Narr))


# random points to evaluate the solution Nmc times for the Monte Carlo quadrature of the error
rps = rand(Nmc).*2.0

# Looping to calcuate the solution and error for all numbers of cells in the list
for k=1:length(Narr)
        println("Solving for k = ", k)
	grid, sol = NTDGD(EulerDGDpe!, Narr[k], tend, uinit, cfl=0.1/(order^2+1), periodic=:true, hot=:true)
	gridus, solus = NTDGD(EulerDG!, Narr[k], tend, uinit, cfl=0.1/(order^2 +1), periodic=:true, hot=:true)
	# Calculation of the error using MC integration
	for i = 1:Nmc
                errarr[k] = errarr[k] + sum(abs.(SolEval(grid, sol(tend), rps[i]) - refsol(rps[i], tend)))
                errarr2[k] = errarr2[k] + sum((SolEval(grid, sol(tend), rps[i]) - refsol(rps[i], tend)).^2)
		# for reference purposes: the unmodified solver
		errarrUS[k] = errarrUS[k] + sum(abs.(SolEval(gridus, solus(tend), rps[i]) - refsol(rps[i], tend)))
                errarr2US[k] = errarr2US[k] + sum((SolEval(gridus, solus(tend), rps[i]) - refsol(rps[i], tend)).^2)
	end
	errarr[k] = errarr[k] / Nmc
	errarr2[k] =sqrt(errarr2[k] / Nmc)
	errarrUS[k] = errarrUS[k] / Nmc
        errarr2US[k] =sqrt(errarr2US[k] / Nmc)
	GC.gc()
end

# Set increased font size for readability
fontsize_theme = Theme(fontsize = 24, colormap=:grayC)
set_theme!(fontsize_theme)
# set axis labels
axis = (xlabel = "Cells", ylabel = "Error", yscale=log10, xscale=log2)
CairoMakie.activate!()

# plot error
fig = scatter(Narr, errarr, axis = axis, color=:black, label="DDG")
scatter!(Narr, errarrUS, color=:black, label="DG", marker=:cross)
# and lines with corresponding orders
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])], label = "order 1, 5", color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^2], label = "order 2, 6", color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^3], label = "order 3, 7", color=:black, linestyle=:dash)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^4], label = "order 4, 8", color=:black)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^5], color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^6], color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^7], color=:black, linestyle=:dash)


lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^8], color=:black)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^9], color=:black, linestyle=:dot)
if order > 7
	lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^10], label = "order 10", color=:black, linestyle=:dashdot)
	lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^11], label = "order 11", color=:black, linestyle=:dash)
	lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^12], label = "order 12", color=:black)
end




# add legend
axislegend(position=:lb)
save("figures/convana1p"*string(order)*"t"*string(tend)*".pdf", fig)

fig = scatter(Narr, errarr2, axis=axis, color=:black, label="DDG")
scatter!(Narr, errarr2US, color=:black, label= "DG", marker=:cross)

lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])], label = "order 1, 5", color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^2], label = "order 2, 6", color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^3], label = "order 3, 7", color=:black, linestyle=:dash)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^4], label = "order 4, 8", color=:black)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^5], color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^6], color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^7], color=:black, linestyle=:dash)


lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^8], color=:black)

	lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^9], color=:black, linestyle=:dot)
if order > 7
	lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^10], label = "order 10", color=:black, linestyle=:dashdot)
	lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^11], label = "order 11", color=:black, linestyle=:dash)
	lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^12], label = "order 12", color=:black)

end

axislegend(position=:lb)
save("figures/convana2p"*string(order)*"t"*string(tend)*".pdf", fig)

GC.gc()
