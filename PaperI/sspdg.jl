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

# Initial conditions for later tests
u1f(x) = sin(pi*x) + 0.5
u2f(x) = x < 1.0 ? -x : 1.0 - (x - 1)
ctestamp = 1/50
u3f(x) = 1 + sin(pi*x)*ctestamp

# flux, entropy flux and entropy
f(u) = u^2/2
F(u) = u^3/3
U(u) = u^2/2

# the flux and entropy derivatives and the entropy Hessian.
df(u) = u
dU(u) = u
ddU(u) = 1

# uncomment to test the quartic entropy.
#U(u) = u^4/4
#F(u) = u^5/5
#dU(u) = u^3
#ddU(u) = 3*u^2

# A Lipschitz-bound for the derivative of the entropy.
Lu(u) = maximum(abs.(ddU.(u)))

# the used polynomial order
order =  6

# the polynomial order identifiable using the collocation points of the error estimate
Horder = order + 1

# collocation nodes for the DG method
nodes, weights = gausslobatto(order+1)

# collocation nodes for the error estimate, when Gauß-Lobatto nodes are used
nodesHO, wHO = gausslobatto(Horder+1)

# collocation nodes for the error estimate, when Gauß-Legendere nodes are used
inodes, weights = gausslegendre(order+1)
inodesHO, wHO = gausslegendre(Horder+1)



x = nodes
xHO = nodesHO
# Matrix evaluates a Ansatz in the modal basis at the nodes
Pm = PaMat(x)
PmHO = PaMat(xHO)

# Modal mass matrix
Mm = Mat(Me, order+1)
MmHO = Mat(Me, Horder+1)

# Modal differentiation matrix
Sm = Mat(Se, order+1)
SmHO = Mat(Se, Horder+1)

# Matix transforms from nodal presentation to modal representation
Om = inv(Pm)
OmHO = inv(PmHO)

# Calculation of the nodal representation of all needed matrices
Mref = transpose(Om)*Mm*Om
Sref = transpose(Om)*Sm*Om

# High order versions
MrefHO = transpose(OmHO)*MmHO*OmHO
SrefHO = transpose(OmHO)*SmHO*OmHO

# Interpolation Matrices, i.e. that embedd the solution into 
# a nodal solution one polynomial degree higher
ToHO = PmHO*I[1:Horder+1, 1:order+1]*Om
ToLO = Pm*I[1:order+1, 1:Horder+1]*OmHO

# Matrices that evaluate the Ansatz at the embedded Gauß-Legendere Nodes
Im = PaMat(inodes)
ImHO = PaMat(inodesHO)

# Matrices that evaluate the Ansatz at the interrior nodes for the error estimate
ToI = Im*Om
ToIHO = ImHO*OmHO

# Matrices that transform from the inner nodal representation to the nodal basis
ToN = inv(ToI)
ToNHO = inv(ToIHO)


# Lax-Friedrichs flux with viscosity mu
function hBurgersLF(ul, ur, mu)
        return 0.5*((ul*ul+ur*ur)/2.0 - mu*(ur-ul))
end

# Local Lax-Friedrichs flux
function hBurgersLLF(ul, ur)
        mu = max(abs(ul), abs(ur))
        return hBurgersLF(ul,ur, mu)
end

# Modified Lax-Friedrchs flux
function hBurgersMLF(ul, ur, mu)
	meanu = 0.5*(ul + ur)
	return meanu.^2 ./2.0 - mu*(ur - ul)./2.0
end
function hBurgersLMLF(ul, ur)
	mu = max(abs(ul), abs(ur))
        return hBurgersMLF(ul,ur, mu)
end


# Numerical entropy flux for the Lax-Friedrichs flux
function HBurgersLF(ul, ur, mu)
	return 0.5*(F(ul) + F(ur) - mu*(U(ur) - U(ul)))
end

# Numerical entropy flux for the local Lax-Friedrichs flux
function HBurgersLLF(ul, ur)
	mu = max(abs(ul), abs(ur))
	return HBurgersLF(ul, ur, mu)
end

# Numerical solution to the one-dimensional Riemann problem
# for a scalar conservation law with conves flux
# xdbt is the position in the fan, i.e. x divided by t
# the a and a inv shall be the derivative of the flux and
# its inverse function
function Riemann1D(xdbt, ul, ur, a, ainv, f)
        s = (f(ul) - f(ur))/(ul-ur)
        if ul > ur
                if xdbt < s
                        return ul
                else
                        return ur
                end
        else
                if xdbt <= a(ul)
                        return ul
                elseif xdbt >= a(ur)
                        return ur
                else
                        return ainv(xdbt)
                end
        end
end

# Solution to the Riemann problem for Burgers' equation
function RiemannBurg(xdbt, ul, ur)
        f(u) = u^2/2
        a(u) = u
        ainv(y) = y
        return Riemann1D(xdbt, ul, ur, a, ainv, f)
end

# Godunov flux for Burgers' equation
function hBurgersGod(ul, ur)
        f(u) = u^2/2
        return f(RiemannBurg(0, ul, ur))
end

# Implementation of a Vanilla DG method for 
# Burgers equation
function BurgersDG!(du, u, p, t)
	dx = p[1]
        dt = p[2]
        Mi = p[3]
        S = p[4]
	N, K = size(u)
        flux(ul, ur) = hBurgersLLF(ul, ur)
	for k=1:N
                v = S*f.(u[k, :])
                v[1] = v[1] + flux(u[mod1(k-1, N), K], u[k, 1])
                v[K] = v[K]  - flux(u[k, K], u[mod1(k+1, N), 1])
		du[k, :] = Mi*v
	end
end

# Helper function for the discrete limiter function that calculates the error bounds
# for the ansatz vector u in nodal representation
function SE(u, p)
        dx = p[1]
        dt = p[2]
        Mi = p[3]
        S = p[4]
        MiHO = p[5]
        SHO = p[6]
        M = p[7]
        MHO = p[8]
        flux(ul, ur) = hBurgersLLF(ul, ur)
        N, K = size(u)
	vlf = zeros(K)
        delta = zeros(N)
        for k=1:N
                # Low order
                v = S*f.(u[k, :])
                v[1] = v[1] + flux(u[mod1(k-1, N), K], u[k, 1])
                v[K] = v[K]  - flux(u[k, K], u[mod1(k+1, N), 1])
                du = Mi*v
                um = ToHO*u[k, :]
                dusHO = df.(ToIHO*um).*ToIHO*MiHO*transpose(SHO)*um
		vlf[1] = -(f(u[k, 1]) - flux(u[mod1(k-1, N), K], u[k, 1]))
                vlf[K] = f(u[k, K]) - flux(u[k, K], u[mod1(k+1, N), 1])
                dusHO = dusHO + ToIHO*ToHO*vlf
                diffHO = dusHO + ToIHO*ToHO*du
                delta[k] = sqrt(dot(wHO, diffHO.^2))
        end
        return delta
end

# Limiter function implementing the approximate solution of the
# stated optimization problem after a full RK step using a 
# (severely) early stopped gradient descent
function DDGlimiters!(u, integ, p, t)
        dx = p[1]
        dt = p[2]
        Mi = p[3]
        S = p[4]
        MiHO = p[5]
        SHO = p[6]
        M = p[7]
        MHO = p[8]
        N, K = size(u)
        du = zeros(N, K)
        du1 = zeros(N,K)
        du2 = zeros(N,K)
        sol = integ.sol
        nit = 5
        
        up = integ.uprev
        integ.f(du, up, p, t) 
        u1 = up + dt*du
        integ.f(du1, u1, p, t)
        u2 = up + 0.25*dt*du + 0.25*dt*du1
        integ.f(du2, u2, p, t)
        delta0 = SE(up, p)
        delta1 = SE(u1, p)
        delta2 = SE(u2, p)

        for k=1:N
                delta = dt*(delta0[k] + delta1[k] + 4*delta2[k])/6.0

		# calculation of the roundoff error
		leps = eps(norm(du[k, :]))
		# as in the semi-discrete case (look below), the correction is only made if
		# it is smaller than the square root of the rounding error, as the calculation 
		# of the error estimate is only well posed if the error is in the magnitude of
		# the square root of the rounding error, or above
		if delta > dt*sqrt(leps)
			for it=1:nit
                		gradU = -dU.(u[k, :])
				gradU = gradU - dot(ones(order+1),Mref*gradU)*ones(order+1)/(dot(ones(order+1), Mref*ones(order+1)))
	                	ngradU = gradU./(sqrt(abs(transpose(gradU)*Mref*gradU)) + 10E-30) 
				unc = delta 
				mdiff = minimum(abs.((Om*u[k, :])[1] .- u[k, :])./(abs.(ngradU) .+ 10E-30))
				eps = min(unc/nit, mdiff)
				u[k, :] = u[k, :] .+ eps*ngradU
			end
		end
        end
end


# semidiscrete DG discretisation,
# stabilizied using Dafermos' entropy rate criterion
# augmented via an error estimate
function BurgersDGDc!(du, u, p, t)
        dx = p[1]
        dt = p[2]
        Mi = p[3]
        S = p[4]
        MiHO = p[5]
        SHO = p[6]
        fac = p[7]
        M = Mref.*dx./2.0

        N, K = size(u)
	# used numerical flux
       	flux(ul, ur) = hBurgersLLF(ul, ur)
	vlf = zeros(K)
	
        for k=1:N
		# calculation ofthe vaniall DG scheme
                v = S*f.(u[k, :])
                um = ToHO*u[k, :]
		v[1] = v[1] + flux(u[mod1(k-1, N), K], u[k, 1])
		v[K] = v[K]  - flux(u[k, K], u[mod1(k+1, N), 1])
                du[k, :] = Mi*v

		# calculation f the error estimate
		ui = ToIHO * um
                dusHO = df.(ui).*ToIHO*MiHO*transpose(SHO)*um
		vlf[1] = -(f(u[k, 1]) - flux(u[mod1(k-1, N), K], u[k, 1]))
		vlf[K] = f(u[k, K]) - flux(u[k, K], u[mod1(k+1, N), 1])
		dusHO = dusHO + ToIHO*ToHO*vlf
                diffHO = dusHO + ToIHO*ToHO*du[k, :]

                diffU = ToHO*dU.(u[k, :]) - dU.(um)
  
		# L2 Norm of the error
  		delta = sqrt(dot(wHO, diffHO.^2))

		# Error in U
                deltaU = sqrt(abs(transpose(diffU)*MrefHO*diffU))

		# calculation of the descent direction
		gradU = - dU.(u[k, :])
                mfgradU = gradU - dot(ones(order+1),Mref*gradU)*ones(order+1)/(dot(ones(order+1), Mref*ones(order+1)))
            
		# Normin of the descent direction
		ngradU = mfgradU./sqrt(abs(transpose(mfgradU)*Mref*mfgradU)+1.0E-30)
                unc = delta*transpose(mfgradU)*M*mfgradU/(dot(mfgradU, M*gradU))
		# calculting upper limit for a sign change
                mdiff = minimum(abs.(((Om*u[k, :])[1] .- u[k, :])./(abs.(ngradU) .+ 1.0E-30)))
		h = min(unc, mdiff/dt)

		# Calculating the size of roundoff errors
		leps = eps(norm(du[k, :]))
		
		# Application of the correction. The correction is only applied if
		# the error estimate is bigger than the square root of the rounding errors
		# because the error estimate itself can be only calculated well posed up to
		# the square root of the rounding error
		if h > sqrt(leps)
                	du[k, :] = du[k, :] + fac*h*ngradU
		end
        end
end


# This function calculates the violation of the 
# semidiscrete entropy inequality at time t
# using the factor $fac$ the error estimate can be scaled
function BurgersdE(grid, sol, t, fac, cfl)
	K, L = size(sol(0.0))
	dx = grid[1, end] - grid[1, 1]
        dt = dx*cfl

	du = zeros(K, L)
        dE = zeros(K)
	
	dF = zeros(K)
	E = zeros(K)
	M = Mref.*dx./2.0
        Mi = inv(Mref).*2.0./dx
        Sl = transpose(Sref)
        MiHO = inv(MrefHO).*2.0./dx
        MHO = MrefHO.*dx./2.0
        SlHO = transpose(SrefHO)
        p = [dx, dt, Mi, Sl, MiHO, SlHO, fac]
        BurgersDGDc!(du, sol(t), p, t)
	for k =1:K
		E1 = dot(ones(order+1), M*U.(sol(t)[k, :]))
	
		u1 = dot(ones(order+1), M*sol(t)[k, :])/dx
		u2 = dot(ones(order+1), M*sol(t+dt)[k, :])/dx
		U1, U2 = U(u1)*dx, U(u2)*dx
		E[k] = E1
		dE[k] = dot(ones(order+1), M*(dU(sol(t)[k, :]).*du[k, :]))#(E2 - E1)/dt 
		dF[k] = (HBurgersLLF(sol(t)[mod1(k-1, K), end],sol(t)[k, 1]) - HBurgersLLF(sol(t)[k, end], sol(t)[mod1(k+1, K), 1]))
	end
	return dE, dF, E
end

# Stub limiter function without effect
function Id!(u,i, p, t)

end

# Driver routine for the time discrete Dafermos Entropy rate
# criterion  stabilized 
# Runge-Kutta discontinuous Galerkin method
# K: Number of Cells, T: Maximum solution time
# u0 initial condition
function NTBurgDiscDG(K, T, u0, cfl = 0.02)
	tspan = (0, T)
        dx = 2/(K)
        dt = dx*cfl
        N = ceil(Int, T/dt)
        grid = dx*collect(1:K) .- dx/2
        gridar = zeros(K, order+1)
        uar = zeros(K, order+1)
        for k=1:order+1
                uar[:, k] = u0.((grid.+ dx/2*x[k]))
                gridar[:, k] = grid.+ dx/2*x[k]
        end
        Mi = inv(Mref).*2.0./dx
        M = Mref.*dx./2.0
        Sl = transpose(Sref)
        MiHO = inv(MrefHO).*2.0./dx
        MHO = MrefHO.*dx./2.0
        SlHO = transpose(SrefHO)
	p = [dx, dt, Mi, Sl, MiHO, SlHO, M, MHO]
        prob = ODEProblem(BurgersDG!, uar, tspan, p)

        sol = solve(prob, SSPRK33(Id!,  DDGlimiters!), dt=dt)
        gridar, sol
end


# Driver routine for numerical tests using the semidiscrete version
# of the Dafermos Entropy rate criterion stabilized DG method
# K is the number of Cells, T the maximum running time, u0 the
# initial condition, fac allows to scale the error estimate
function NTBurgDGSSP(K,T, u0, fac = 1.0, cfl = 0.02)
        # initialization of auxillary variables
	tspan = (0, T)
        dx = 2/(K)
        dt = dx*cfl # assuming cmax = 1
	N = ceil(Int, T/dt)
        grid = dx*collect(1:K) .- dx/2
        gridar = zeros(K, order+1)
	uar = zeros(K, order+1)
        # especially the grid and initial condition
	for k=1:order+1
                uar[:, k] = u0.((grid.+ dx/2*x[k]))
                gridar[:, k] = grid.+ dx/2*x[k]
        end
        
        # Scaling of the mass and differentiation matrices
	Mi = inv(Mref).*2.0./dx
        M = Mref.*dx./2.0
        Sl = transpose(Sref)
        MiHO = inv(MrefHO).*2.0./dx
        MHO = MrefHO.*dx./2.0
        SlHO = transpose(SrefHO)
	
        # vector of paramters
	p = [dx, dt, Mi, Sl, MiHO, SlHO, fac]
	prob = ODEProblem(BurgersDGDc!, uar, tspan, p)
	sol = solve(prob, SSPRK33(), dt=dt, adaptive=false)
	gridar, sol
end


# The following two functions plot a solution sol on a grid g at time t
# the title in the legend is name, the solution can be overlayed onto a 
# previus one if overlay is specified. Linestyle and some space around, as bd, can be specified
# the second version of the function is needed to provide animations
function PlotSol(g, sol, t; name = "solution", overlay = false, ls = :solid, bd = 0.3)
	N, K = size(sol(0.0)[:, :])
	upper = bd + maximum(sol(0.0))
	lower = minimum(sol(0.0)) - bd
        axis = (xlabel = "t", ylabel = "u(x, t)")
	if overlay
		p = lines!(g[1, :], sol(t)[1, :], color=:black, linestyle = ls, label = name)
	else
	        p = lines(g[1, :], sol(t)[1, :], axis=axis, color=:black, linestyle = ls, label = name)
        end

        for k=2:N
		lines!(g[k, 1:end], sol(t)[k, 1:end], color=:black, linestyle = ls)
	end
        #axislegend()	
        return p
end

function PlotOSol(g, sol, to; name = "solution", overlay = false, ls = :solid, bd = 0.3)
        N, K = size(sol(0.0)[:, :])
        upper = bd + maximum(sol(0.0))
        lower = minimum(sol(0.0)) - bd
        axis = (xlabel = "t", ylabel = "u(x, t)")
        if overlay
		yv = @lift sol($to)[1, :]
                p = lines!(g[1, :], yv, color=:black, linestyle = ls, label = name)
        else
		yv = @lift sol($to)[1, :]
                p = lines(g[1, :], yv, axis=axis, color=:black, linestyle = ls, label = name)
        end

        for k=2:N
		yv = @lift sol($to)[k, :]
                lines!(g[k, 1:end], yv, color=:black, linestyle = ls)
        end
        axislegend()
        return p
end

# Function saves an animation of the solution sol on grid g under the filename
# fname. The solution is evaluated up to tmax and using steps intermediate evaluations
function MKanim(g, sol, steps, tmax, fname)
	nframes = steps
	framerate = 30
	to = Observable(0.0)
	tit = range(0.0, tmax, length=steps)
	fig, ax, lineplot = PlotOSol(g, sol, to)
	record(fig, fname, tit; framerate=framerate) do t
		to[]=t
	end
end
