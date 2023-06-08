# File performs a convergence analysis for the discrete scheme
# The initial condition sin(\pi x)/50 + 1 is used
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
	return SCeval(grid[lb, :], sol[lb, :], x)
end

# solving for the initial condition
# uinit(x) = sin(pi*x) / 50 .+ 1.0
function refsol(x, t, N = 100)
        ures = 1.0
        for k = 1:N
                xpred = x - t*ures
                ures = u3f(xpred)
        end
        return ures
end

# load the solver
include("sspdg.jl")

# Number of random evaluation points
# where the norm of the error is calcualted
Nmc = 10000

# end time where the solution error is calculated
tend = 8.0

# Array of possible cell numbers
Narr = [10, 15, 20, 25, 30, 40, 50]
cflp = 0.02*Narr[1]
# Arrays for saving the errors
errarr = zeros(length(Narr))
errarr2 = zeros(length(Narr))

# random points for MC integration of the error
rps = rand(Nmc).*2.0


for k=1:length(Narr)
	grid, sol = NTBurgDiscDG(Narr[k], tend, u3f, cflp/Narr[k])
	for i = 1:Nmc
		errarr[k] = errarr[k] + abs(SolEval(grid, sol(tend), rps[i]) - refsol(rps[i], tend))
		errarr2[k] = errarr2[k] + (SolEval(grid, sol(tend), rps[i]) - refsol(rps[i], tend))^2

	end
	errarr[k] = errarr[k] / Nmc
	errarr2[k] =sqrt(errarr2[k] / Nmc)

end

fontsize_theme = Theme(fontsize = 28, colormap=:grayC)
set_theme!(fontsize_theme)
axis = (xlabel = "Cells", ylabel = "Error", yscale=log10, xscale=log10)
CairoMakie.activate!()
fig = scatter(Narr, errarr, axis = axis, color=:black)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])], label = "order 1", color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^2], label = "order 2", color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^3], label = "order 3", color=:black, linestyle=:dash)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^4], label = "order 4", color=:black)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^5], label = "order 5", color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^6], label = "order 6", color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr[1], errarr[1]*(Narr[1]/Narr[end])^7], label = "order 7", color=:black, linestyle=:dash)

axislegend(position=:lb)
save("dconvana1p"*string(order)*"amp"*string(Int(100*ctestamp))*"t"*string(Int(tend))*".pdf", fig)

fig = scatter(Narr, errarr2, axis=axis, color=:black)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])], label = "order 1", color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^2], label = "order 2", color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^3], label = "order 3", color=:black, linestyle=:dash)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^4], label = "order 4", color=:black)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^5], label = "order 5", color=:black, linestyle=:dot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^6], label = "order 6", color=:black, linestyle=:dashdot)
lines!([Narr[1], Narr[end]], [errarr2[1], errarr2[1]*(Narr[1]/Narr[end])^7], label = "order 7", color=:black, linestyle=:dash)

axislegend(position=:lb)
save("dconvana2p"*string(order)*"amp"*string(Int(100*ctestamp))*"t"*string(Int(tend))*".pdf", fig)


