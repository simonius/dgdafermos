# inclusion of needed packages
# Solving the ODE systems
using DifferentialEquations
# Inverting the needed matrices
using LinearAlgebra
# Plotting
using Makie
using CairoMakie
using GLMakie
using LaTeXStrings
# Generation of Gauss quadrature nodes
using FastGaussQuadrature

# The file bases contains functions that implement the matrix elements
# of the needed matrices
include("bases.jl")

# Dissipation Operators
include("DispMatGen.jl")

# flux of Burgers' equation
f(u) = u.^2 ./2
# flux derivative
df(u) = u
# Entropy
U(u) = u[1]^2
# entropy variables
dU(u) = 2*u
# entropy flux
F(u) = 2*u[1]^3 ./3

# number of conserved quantities
ncons = 1

# initial conditions
u1f(x) = [sinpi(x)+0.5]
u2f(x) = [x < 1.0 ? -x : 2-x]
u3f(x) = [x < 1.3 ? 0.0 : 1.0]

# Speed estimate for HLL type numerical flux.
# eps is used to remove a divide by zero.
function HLLfopnorm(ul, ur)
	al = min(ul[1], ur[1])
        ar = max(ul[1], ur[1])
	eps = 1.0E-6
	return al - eps, ar + eps
end

# simplified estimate for LLF based fluxes
function LLFfopnorm(ul, ur)
	c = max(abs(ul[1]), abs(ur[1]))
	return c
end

# HLL numerical flux function
function hBurgersHLL(ul, ur)
	al, ar = HLLfopnorm(ul, ur)
	if (0 < al)
                return f(ul)
        elseif (0 < ar)
		return (al*ar*(ur-ul) + ar*f(ul) - al*f(ur)) / (ar-al)
        else
		return f(ur)
        end
end

# HLL numerical entropy flux
function HBurgersHLL(ul, ur)
	al, ar = HLLfopnorm(ul, ur)
	
	um = um = (ar*ur - al*ul)/(ar-al) + (f(ul) - f(ur))/(ar-al)
        if (0 < al)
		flb = al*U(ul) + (ar - al) * U(um) - ar*U(ur) + F(ur)
                fub = F(ul)
                Hhll = fub
	elseif (0 < ar)
		flb = ar*U(um) - ar*U(ur) + F(ur)
                fub = al*U(um) - al*U(ul) + F(ul)
		Hhll = 0.5*(flb + fub)
	else
		flb = F(ur)
                fub = F(ul) - (ar - al)*U(um)+ ar*U(ur) - al*U(ul)
                Hhll = flb
	end
	return Hhll
end

# HLL type entropy inequality predictor
function HLLDispOp(ul, ur)
	al, ar = HLLfopnorm(ul, ur)
        
	um = (ar*ur - al*ul)/(ar-al) + (f(ul) - f(ur))/(ar-al)
        lUd = (ar - al)*U(um) + al*U(ul) - ar*U(ur) - (F(ul) - F(ur))
        # We do not know in which cell the dissipation takes place
	return lUd
end

# local Lax-Friedrichs flux
function hBurgersLLF(ul, ur)
	c = LLFfopnorm(ul, ur)
	return 0.5*(f(ul) + f(ur) + c*(ul - ur))
end

# local Lax-Friedrichs entropy flux
function HBurgersLLF(ul, ur)
	c = LLFfopnorm(ul, ur)
	return 0.5*(F(ul) + F(ur) + c*(U(ul) - U(ur)))
end

# local Lax-Friedrichs based entropy inequality predictor
function LLFDispOp(ul, ur)
	c = LLFfopnorm(ul, ur)
	um = 0.5*(ul + ur) + (f(ul) - f(ur))/(2*c)
        return c*(2*U(um) - U(ul) - U(ur)) - F(ul) + F(ur)
end

# Selection of the numerical flux, numerical entropy flux and entropy inequality predictor
# possible choices: hBurgersHLL, hBurgersLLF
flux = hBurgersHLL
# possible choices: HBurgersHLL, HBurgersLLF
Flux = HBurgersHLL
# possible choices: HLLDispOp,  LLFDispOp 
dispOp = HLLDispOp

# the simulation domain
xr = 2.0
xl = 0.0

# the used polynomial order
order = 6

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

# Projects an initial condition onto the Legendre polynomials
# used to discretise the initial condition
function L2Proj(uf, Nmax)
        c = zeros(Nmax)
        for k=1:Nmax
                integ(x) = uf(x)*Pleg(k-1, x)
                int, err = quadgk(integ, -1.0, 1.0, maxevals=100)
                c[k] = (2*k-1)*int*0.5
        end
        return c
end


# DG discretisation of the
# Burgers equation
# u is the N x K x 1 tensor of conserved variables,
# for N cells, K nodes and 1 conserved variable
function BurgersDG!(du, u, p, t)
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

# A maximal dissipation estimate, hardened using the projection onto lower order
# polynomials in the case order > 2.
# Per Edge, i.e. dE refers to the dissipation taking place at the cell boundary k+1
function DispEstPE(u)
        N, K, L = size(u)
        dE = zeros(N)
        for n=1:N
                ul = u[n, :, :]
                ur = u[mod1(n+1, N), :, :]
                dEloc = dispOp(ul[K, :], ur[1, :])
                dE[n] = dEloc
        end


        if order > 2 
                ord = K-1
                for n=1:N
                	uLOl = ToLOset[ord]*u[n, :, :]
                        uLOr = ToLOset[ord]*u[mod1(n+1, N), :, :]
                        dEloc = dispOp(uLOl[end, :], uLOr[1, :])
                        dE[n] = min(dE[n], dEloc)
                end
        end
        return dE
end

# The total entropy functional
function Efunc(u, par)
        dx = par[1]
        K, L = size(u)
        Uv = zeros(K)
        for k=1:K
                Uv[k] = U(u[k, :])
        end
        return dx*0.5*dot(weights, Uv)
end

# A DG discretisation of the Euler equations, corrected in entropy by a maximal
# entropy rate estimate per cell edge.
function BurgersDGDpe!(du, u, par, t)
        N, K, L = size(u)
        dx, dt = par[1], par[2]
        BurgersDG!(du, u, par, t)

        dEPE = DispEstPE(u)
        dEDG = zeros(N)
        dEDGpe = zeros(N)
        dEdd = zeros(N)
        dd = zeros(N, K, L)
        dUdtDG = zeros(K)
        dUdtdd = zeros(K)
        alp = zeros(N)
        leps = sqrt(eps(1.0))
        for n=1:N
                dd[n, :, :] = G*u[n, :, :]
                for k=1:K
                        dUdtDG[k] = dot(dU(u[n, k, :]), du[n, k, :])
                        dUdtdd[k] = dot(dU(u[n, k, :]), dd[n, k, :])
                end
                dEDG[n] = dot(weights*dx/2.0, dUdtDG)
                dEdd[n] = Efunc(u[n, :, :] + G*u[n, :, :], par) - Efunc(u[n, :, :], par)
                dEDG[n] = dEDG[n] - (Flux(u[mod1(n-1, N), K, :], u[n, 1, :]) - Flux(u[n, K, :], u[mod1(n+1, N), 1, :]))

                if dEdd[n] > leps
                        println("antidissipative with ", dEdd[n], " ", leps2)
                end
                # the needed alpha for entropy stability
                alpES = max(0.0, -dEDG[n]*dEdd[n]/(dEdd[n]^2 + leps^2))
                alp[n] = alp[n] + alpES
                # update the entropy error of the scheme - now it will be non-positive.
                dEDG[n] = dEDG[n] + alpES*dEdd[n]
                # splitting of the dissipation per edge.
                dEDGl = dEDG[n]*(dEPE[mod1(n-1, N)]/(dEPE[mod1(n-1, N)] + dEPE[n] - leps))
                dEDGr = dEDG[n]*(dEPE[n]/(dEPE[mod1(n-1, N)] + dEPE[n] - leps))
                dEDGpe[n] = dEDGpe[n] + dEDGr
                dEDGpe[mod1(n-1, N)] = dEDGpe[mod1(n-1, N)] + dEDGl
        end
        for n=1:N
                # the per unit correction Entropy rate at the edge
                dEeg = dEdd[n] + dEdd[mod1(n+1, N)]
                st = dEPE[n] - 1.0*dEDGpe[n]
                alpER = max(0.0, st*dEeg/(dEeg^2 + leps^2))
                alp[n] = alp[n] + alpER
                alp[mod1(n+1, N)] = alp[mod1(n+1, N)] + alpER
        end
        maxval = 1.0/dt
        for n=1:N
                du[n, :, :] = du[n, :, :] + min(maxval, 1.0*alp[n])*dd[n, :, :]
        end

        if par[end] == :false
                du[end, :, :] .= 0.0
                du[1, :, :] .= 0.0
                alp[1], alp[end] = 0.0, 0.0
        end
        return alp, dEDG, dEPE, dEdd
end

# Driver routine for numerical tests using the semidiscrete version
# of the Dafermos Entropy rate criterion stabilized DG method
# K is the number of Cells, T the maximum running time, u0 the
# initial condition, periodic switches on periodic BCs,
# hot activates a High-Order time integration method for convergence tests.
function NTDGD(solver, K,T, u0; cfl = 0.02, periodic=:true, hot=:false)
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
                for l=1:ncons
                        c = L2Proj(x->u0(grid[n] + dx/2*x)[l], order+1)
                        uar[n, :, l] = Pm*c
                end
                for k=1:order+1
                   #     uar[n, k, :] = u0(grid[n]+ dx/2*x[k])
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
function PlotSol(g, sol, t, comp; name = "solution", overlay = false, ls = :solid, bd = 0.3, yl = "u(x, t)")
        N, K, L= size(sol(0.0)[:, :, :])
        upper = bd + maximum(sol(0.0)[:, :, comp])
        lower = minimum(sol(0.0)[:, :, comp]) - bd
        axis = (xlabel = L"x", ylabel = yl)
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

