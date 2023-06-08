# inclusion of needed packages
# Solving the ODE systems
using DifferentialEquations
# Inverting the needed matrices
using LinearAlgebra
# Plotting
using Makie
using CairoMakie
using GLMakie
# Generation of Gauss quadrature nodes
using FastGaussQuadrature

# The file bases contains functions that implement the matrix elements
# of the needed matrices
include("bases.jl")

# Dissipation Operators
include("DispMatGen.jl")

# auxialry definitions for the euler equations
include("euler.jl")
f = fEuler
df = dfEuler
dU = dUEuler
Ufunc = PUEuler
flux = hEulerLLF
Flux = HEulerLLF
F = PFEuler

xr = 10.0
xl = 0.0
# the used polynomial order
order =  3

# collocation nodes for the DG method
nodes, weights = gausslobatto(order+1)
x = nodes



# Matrix evaluates a Ansatz in the modal basis at the nodes
Pm = PaMat(x)
# Modal mass matrix
Mm = Mat(Me, order+1)
# Modal differentiation matrix
Sm = Mat(Se, order+1)
# Matix transforms from nodal presentation to modal representation
Om = inv(Pm)
# Calculation of the nodal representation of all needed matrices
Mref = transpose(Om)*Mm*Om
Sref = transpose(Om)*Sm*Om
weights = ones(1, order + 1)*Mref

G = DispMats(x, weights)

# Matrices for projection onto a lower order solution
ToLOset = Array{Any}(undef, order)
for ord = 2:order
        nLO, wLO = gausslobatto(ord)
        PmLO = PaMat(nLO)
        ToLO = PmLO*I[1:ord, 1:order+1]*Om
        ToLOset[ord] = ToLO
end

khmax = floor(Int, (order+1)/2)
IEM = zeros(2, order+1)
for k=1:khmax
        IEM[1, k] = Plag(nodes[1:khmax], k, 0.0)
        IEM[2, khmax+k] = Plag(nodes[khmax+1:end], k, 0.0)
end

R = zeros(2, order+1)
for k=1:order
        R[1, k+1] = Plag(nodes[2:end], k, -1.0)
        R[2, k] = Plag(nodes[1:end-1], k, 1.0)
end

# DG discretisation of the
# Euler equations
# u is the N x K x L tensor of conserved variables,
# for N cells, K nodes and L conserved variables
function EulerDG!(du, u, p, t)
        dx = p[1]
        dt = p[2]
        Mi = p[3]
        S = p[4]
        N, K, L = size(u)
        fv = zeros(K, L)
        v = zeros(K, L)
        for n=1:N
                for k=1:K
                        fv[k, :] = f(u[n, k, :])
                end
                for l=1:L
                        v[:, l] = S*fv[:, l]
                end

                v[1, :] = v[1, :] + flux(u[mod1(n-1, N), K, :], u[n, 1, :])
                v[K, :] = v[K, :] - flux(u[n, K, :], u[mod1(n+1, N), 1, :])
                for l=1:L
                        du[n, :, l] = Mi*v[:, l]
                end
        end
        if !p[end]
                du[1, :, :] .= 0.0
                du[end, :, :] .= 0.0
        end
end

# The Dissipation Prediction of a local Lax-Friedrichs style entropy inequality predictor,
# i.e. the HLL predictor from the publication with the trivial speed estimate -a_l = c = a_r
function LLFDispOp(ul, ur)
        pl, pr = p(ul), p(ur)
        amax = max(a(ul, pl), a(ur, pr))
        cbound=  1.0*(max(abs(ul[2]/ul[1]), abs(ur[2]/ur[1])) + amax)
        um = 0.5*(ul + ur) + (fEuler(ul) - fEuler(ur))/(2*cbound)
        return cbound*(2*Ufunc(um) - Ufunc(ul) - Ufunc(ur)) - F(ul) + F(ur)
end

# A maximal dissipation estimate, hardened using the projection onto lower order
# polynomials in the case order > 2.
# Per Edge, i.e. dE refers to the dissipation taking place at the cell boundary k+1
function LLFDispEstPE(u)
        N, K, L = size(u)
        dE = zeros(N)
        for n=1:N
                ul = u[n, :, :]
                ur = u[mod1(n+1, N), :, :]
                dEloc = LLFDispOp(ul[K, :], ur[1, :])
                dE[n] = dEloc
        end


        if order > 2
                ord = K-1
                for n=1:N
                        uLOl = ToLOset[ord]*u[n, :, :]
                        uLOr = ToLOset[ord]*u[mod1(n+1, N), :, :]
                        dEloc = LLFDispOp(uLOl[end, :], uLOr[1, :])
                        dE[n] = min(dE[n], dEloc)
                end
        end
        return dE
end

# A DG discretisation of the Euler equations, corrected in entropy by a maximal
# entropy rate estimate per cell edge.
function EulerDGDpe!(du, u, par, t)
        N, K, L = size(u)
        dx, dt = par[1], par[2]
        EulerDG!(du, u, par, t)

        dEPE = LLFDispEstPE(u)
        dEDG = zeros(N)
        dEDGpe = zeros(N)
        dEdd = zeros(N)
        dd = zeros(N, K, L)
        dUdtDG = zeros(K)
        dUdtdd = zeros(K)
        alp = zeros(N)
        leps1 = sqrt(eps(1.0))
        leps2 = 10^4*eps(1.0)
        for n=1:N
                dd[n, :, :] = G*u[n, :, :]
                for k=1:K
                        dUdtDG[k] = dot(dU(u[n, k, :]), du[n, k, :])
                        dUdtdd[k] = dot(dU(u[n, k, :]), dd[n, k, :])
                end
                dEDG[n] = dot(weights*dx/2.0, dUdtDG)
                dEdd[n] = dot(weights*dx/2.0, dUdtdd)
                dEDG[n] = dEDG[n] - (Flux(u[mod1(n-1, N), K, :], u[n, 1, :]) - Flux(u[n, K, :], u[mod1(n+1, N), 1, :]))

                if dEdd[n] > leps2
                        println("antidissipative with ", dEdd[n], " ", leps2)
                        leps2 = 2*(dEdd[n] + leps2)
                end
                # the needed alpha for entropy stability
                alpES = max(0.0, -dEDG[n]*dEdd[n]/(dEdd[n]^2 + leps2^2))
                alp[n] = alp[n] + alpES
                # update the entropy error of the scheme - now it will be non-positive. 
                dEDG[n] = dEDG[n] + alpES*dEdd[n]
 		# splitting of the dissipation per edge.
                dEDGl = dEDG[n]*(dEPE[mod1(n-1, N)]/(dEPE[mod1(n-1, N)] + dEPE[n] - leps1))
                dEDGr = dEDG[n]*(dEPE[n]/(dEPE[mod1(n-1, N)] + dEPE[n] - leps1))
                dEDGpe[n] = dEDGpe[n] + dEDGr
                dEDGpe[mod1(n-1, N)] = dEDGpe[mod1(n-1, N)] + dEDGl
        end
        for n=1:N
                # the per unit correction Entropy rate at the edge
                dEeg = dEdd[n] + dEdd[mod1(n+1, N)]
                st = dEPE[n] - dEDGpe[n]
                alpER = max(0.0, st*dEeg/(dEeg^2 + leps1^2))
                alp[n] = alp[n] + alpER
                alp[mod1(n+1, N)] = alp[mod1(n+1, N)] + alpER
        end
        maxval = 0.5/dt
        for n=1:N
                if alp[n] == maxval
                        println("to high lam needed")
                end
                du[n, :, :] = du[n, :, :] + min(maxval, alp[n])*dd[n, :, :]
        end

        if par[end] == :false
                du[end, :, :] .= 0.0
                du[1, :, :] .= 0.0
                alp[1], alp[end] = 0.0, 0.0
        end
        return alp, dEDG, dEPE, dEDGpe
end




# Driver routine for numerical tests using the semidiscrete version
# of the Dafermos Entropy rate criterion stabilized DG method
# K is the number of Cells, T the maximum running time, u0 the
# initial condition, periodic switches on periodic BCs, 
# hot activates a High-Order time integration method for convergence tests.
function NTDGD(solver, K,T, u0; cfl = 0.02, periodic=:false, hot=:false)
        # initialization of auxillary variables
        tspan = (0, T)
        dx = (xr-xl)/(K)
        dt = dx*cfl # assuming cmax = 1
        N = ceil(Int, T/dt)
        grid = xl .+ dx*collect(1:K) .- dx/2
        gridar = zeros(K, order+1)
        uar = zeros(K, order+1, ncons)
        # especially the grid and initial condition
        for n=1:K
                for k=1:order+1
                        uar[n, k, :] = u0(grid[n]+ dx/2*x[k])
                        gridar[n, k] = grid[n] + dx/2*x[k]
                end
        end

        # Scaling of the mass and differentiation matrices
        Mi = inv(Mref).*2.0./dx
        M = Mref.*dx./2.0
        Sl = transpose(Sref)

        # vector of paramters
        p = [dx, dt, Mi, Sl, M, Mi*Sref, G, periodic]
        prob = ODEProblem(solver, uar, tspan, p)
        if !hot
                sol = solve(prob, SSPRK43(), dt=dt, adaptive=false)
        else
                sol = solve(prob, DP8(), dt=dt, adaptive=false)
        end

        gridar, sol
end

# The following two functions plot a solution sol on a grid g at time t
# the title in the legend is name, the solution can be overlayed onto a
# previus one if overlay is specified. Linestyle and some space around, as bd, can be specified
# the second version of the function is needed to provide animations
function PlotSol(g, sol, t, comp; name = "solution", overlay = false, ls = :solid, bd = 0.3)
        N, K, L= size(sol(0.0)[:, :, :])
        upper = bd + maximum(sol(0.0)[:, :, comp])
        lower = minimum(sol(0.0)[:, :, comp]) - bd
        axis = (xlabel = "x", ylabel = "u(x, t)")
        if overlay
                p = lines!(g[1, :], Float32.(sol(t)[1, :, comp]), color=:black, linestyle = ls, label = name)
        else
                p = lines(g[1, :], Float32.(sol(t)[1, :, comp]), axis=axis, color=:black, linestyle = ls, label = name)
        end

        for k=2:N
                lines!(g[k, 1:end], Float32.(sol(t)[k, 1:end, comp]), color=:black, linestyle = ls)
        end
        return p
end

