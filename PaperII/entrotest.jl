# File carries out an entropy analysis, i.e.
# the classical entropy in-equality and Dafermos'
# principle of maximal entropy dissipation is tested
# for a shock calculation

# load solvers
include("eulerdg.jl")
include("lfsolver.jl")
using CairoMakie

# end time of the calculation
tmax = 2.0
# number of instants in time where the violation of the entropy
# equality is plotted as a heatmap
tsteps = 1000
# Number of cells used for the test (100 in production)
N = 100
# number of cells used in the reference MUSCL solver (30000 in production)
NLF = 30000
# calculate DG solutions that should be analysed
grid, sol =  NTDGD(EulerDGDpe!, N, tmax, shuosh6a, cfl = 0.2/order^2)

# Matrices to save the entropy
Emat = zeros(N, tsteps)
# and the entropy dissipation
dEmat = zeros(N, tsteps)
# vector to save the total entropy
Etot = zeros(tsteps)
# and the derivative of the total entropy
dEtot = zeros(tsteps)

# set the fontsize of the ploting library for readability.
fontsize_theme = Theme(fontsize = 28, colormap=:grayC)
set_theme!(fontsize_theme)

# array of time instants, as a reference for plots
tar = collect(0.0:tmax/tsteps:tmax)[1:tsteps]

# loop for the calculation of the total entropies 
# and the violation of the entropy equality
for i=1:tsteps
	uDG = sol(tar[i])
	N, K, L = size(uDG)
	Uvec = zeros(K)
	dUdt = zeros(K)
	duDG = zeros(N, K, L)
	EulerDGDpe!(duDG, uDG, sol.prob.p, tar[i])
	dx = sol.prob.p[1]

	for n=1:N
		for k=1:K
			Uvec[k] = Ufunc(uDG[n, k, :])
			dUdt[k] = dot(dU(uDG[n, k, :]), duDG[n, k, :])
		end
		Emat[n, i] = dot(weights, Uvec)/2.0*sol.prob.p[1]
		dEdt = dot(weights*dx/2.0, dUdt)
		dEmat[n, i] = dEdt - (Flux(uDG[mod1(n-1, N), K, :], uDG[n, 1, :]) - Flux(uDG[n, K, :], uDG[mod1(n+1, N), 1, :]))
	end
	Etot[i] = sum(Emat[:, i]) + tar[i]*F(uDG[1, 1, :]) - tar[i]*F(uDG[end, end, :])
end

# array of time instants, as a reference for plots
tar = collect(0.0:tmax/tsteps:tmax)[1:tsteps]

# heatmap of the negative logarithm of the violation of the entropy equality
# the small value 1.0E-50 is added for regularization
fig = Makie.heatmap(2:N-1, tar, log10.(max.(eps(1.0), -dEmat[2:N-1, :])), axis = (xlabel = "Cells", ylabel="Time"))
Colorbar(fig.figure[1, 2], fig.plot)
CairoMakie.activate!()
save("figures/entroana-"*string(order)*".pdf", fig)

#heatmap of the positive logarithm of the violation of the entropy equlity
fig = Makie.heatmap(2:N-1, tar, log10.(max.(eps(1.0), dEmat[2:N-1, :])), axis = (xlabel = "Cells", ylabel="Time"))
Colorbar(fig.figure[1, 2], fig.plot)
CairoMakie.activate!()
save("figures/entroana+"*string(order)*".pdf", fig)

# Array for the total entropy of a fine-grained Godunov reference solution
EtotLF = zeros(tsteps)
# calculation of the reference solution
solLF = EulerLF(shuosh6a, tmax, NLF, 2.5)
tsLF = size(solLF)[1]
dxLF = (xr - xl)/NLF
# filling of the array with total entropies
for i=1:tsteps
	lfti = max(1, floor(Int, i/tsteps*tsLF))
	for n = 1:NLF 
		EtotLF[i] = EtotLF[i]+ Ufunc(solLF[lfti, n, :])*dxLF
	end
	EtotLF[i] = EtotLF[i] + tar[i]*F(solLF[lfti, 1, :]) - tar[i]*F(solLF[lfti, end, :])
end


# figure with the total entropies over time
fig = lines(tar, Etot, axis = (xlabel = "time", ylabel = "total entropy"), label = "total entropy, DDG", linestyle=:solid, color=:black)
#lines!(tar, dEtot, label = "total entropy, DRKDG", linestyle = :dash, color=:black)
lines!(tar, EtotLF, label = "LF total entropy", linestyle=:dot, color=:black)
axislegend(position=:lb)
save("figures/totentroana"*string(order)*".pdf", fig)

GLMakie.activate!()

