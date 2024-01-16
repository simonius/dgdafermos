# Tests for the Buckley-Leverett equation

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

# the continuous flux function
const abl = 1.0
f(u) = u.^2.0./(u.^2.0 .+ abl.*(1.0.-u).^2)
df(u) = 2.0*u.^1.0./(u.^2.0 + (-u .+ 1.0).^2.0) + u.^2.0.*(-2.0*u.^1.0 + 2.0*(-u .+ 1.0).^1.0)./(u.^2.0 + (-u .+ 1.0).^2.0).^2

# The entropy function
U(u) = (u[1] .- 0.5).^4
# and the entropy variables
dU(u) = 4*(u .- 0.5).^3
# conserved variables from entropy variables
u(v) = cbrt(v/4.0) + 0.5

# the number of cubature nodes for the entropy flux functions
Nfq = 30
# the cubature nodes and weights
fqg, fqw = gausslobatto(Nfq)
# rescaling of the nodes and weights to lie in [0.0, 1.0]
fqg = 0.5*(fqg .+ 1.0)
fqw = 0.5*fqw

# the numerical entropy flux, integrated using numerical quadrature
# of the derivative of the entropy flux. 
function F(u)
        return u[1]*dot(fqw, df(u[1]*fqg).*dU(u[1]*fqg))
end

# number of conserved variables
ncons = 1

# the initial Riemann data
u1f(x) = [x < 1.0 ? 1.0 : 0.0]

# calculation of the speed estimate for all HLL type 
# estimates a splitting into the increasing and decreasing parts of f(u) is used.
# eps is used to mitigate division by zero in case the left and right speed coincide
function HLLfopnorm(ul, ur)
        minv = min(ul[1], ur[1])
        maxv = max(ul[1], ur[1])
        al, ar = 0.0, 0.0
        eps = 1.0E-6
        if maxv < 0.5   # df is monotone increasing
                al = df(minv) - eps
                ar = df(maxv) + eps
        elseif minv > 0.5 # df is monotone decreasing
                al = df(maxv) - eps
                ar = df(minv) + eps
        else    # df has a maximum at 0.5
                ar = df(0.5) + eps
                al = min(df(minv), df(maxv)) - eps
        end
        return al, ar
end

# A estimate for local LF type fluxes
function LLFfopnorm(ul, ur)
	c = max(abs(df(ul)[1]), abs(df(ur)[1]))
	return c
end

# the HLL flux for the BL equations
function hBLHLL(ul, ur)
	al, ar = HLLfopnorm(ul, ur)
	if (0 < al)
                return f(ul)
        elseif (0 < ar)
		return (al*ar*(ur-ul) + ar*f(ul) - al*f(ur)) / (ar-al)
        else
		return f(ur)
        end
end

# the HLL numerical entropy flux for the BL equations
function HBLHLL(ul, ur)
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

# a HLL based entropy dissipation estimate for the BL equation
function HLLDispOp(ul, ur)
	al, ar = HLLfopnorm(ul, ur)
        
	um = (ar*ur - al*ul)/(ar-al) + (f(ul) - f(ur))/(ar-al)
        lUd = (ar - al)*U(um) + al*U(ul) - ar*U(ur) - (F(ul) - F(ur))
        # We do not know in which cell the dissipation takes place
	return lUd
end

# A local Lax-Friedrichs flux for the Buckley-Leverett equations
function hBLLLF(ul, ur)
	c = LLFfopnorm(ul, ur)
	return 0.5*(f(ul) + f(ur) + c*(ul - ur))
end

# A local Lax-Friedrichs entropy flux for the BL equations
function HBLLLF(ul, ur)
	c = LLFfopnorm(ul, ur)
	return 0.5*(F(ul) + F(ur) + c*(U(ul) - U(ur)))
end

# A LLF based entropy dissipation estimate for the BL equations
function LLFDispOp(ul, ur)
	c = LLFfopnorm(ul, ur)
	um = 0.5*(ul + ur) + (f(ul) - f(ur))/(2*c)
        return c*(2*U(um) - U(ul) - U(ur)) - F(ul) + F(ur)
end

# Selection of the numerical flux. Possible are hBLHLL and hBLLLF
flux = hBLHLL
# Selection of numerical entropy flux. Possible are hBLHLL and hBLLLF
Flux = HBLHLL
# Selection of the entropy inequality predictor - HLLDispOp or LLFDispOp
dispOp = HLLDispOp

# Right and left end of the domain
xr = 2.0
xl = 0.0

# the used polynomial order
order = 7

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

# Projects the function uf onto the Legendre basis 
# used to project the initial condition
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
# Buckley-Leverett equation
# u is the N x K x 1 tensor of conserved variables,
# for N cells, K nodes and 1 conserved variable
function BLDG!(du, u, p, t)
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

# The total entropy functional for a single cell
function Efunc(u, par)
        dx = par[1]
        K, L = size(u)
        Uv = zeros(K)
        for k=1:K
                Uv[k] = U(u[k, :])
        end
        return dx*0.5*dot(weights, Uv)
end

# A DG discretisation of the Buckley-Leverett equation, corrected in entropy by a maximal
# entropy rate estimate per cell edge.
function BLDGDpe!(du, u, par, t)
        N, K, L = size(u)
        dx, dt = par[1], par[2]
        BLDG!(du, u, par, t)

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
                        println("antidissipative with ", dEdd[n], " ", leps)
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
function PlotSol(g, sol, t, comp; name = "solution", overlay = false, ls = :solid, bd = 0.4, yl = "u(x, t)")
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

