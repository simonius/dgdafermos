using DifferentialEquations
using QuadGK

include("bases.jl")

# The C^\infty Friedrichs test function, used for cutting of the viscosity in the cell
function FB(x)
        if abs(x) < 1.0
                y =  exp(-1/(1-x^2))*exp(1.0)
                return y
        else
                return 0.0
        end
end

# Solves the Heat equation using a Galerkin ansatz on the cell
function SolHESGal(u0ar, tmax)
        D = inv(Mref)*Sref
        L = zeros(order+1, order+1)
        NB = diagm(ones(order+1))
        Visc = FB
        for k=1:order+1
                for l=1:order+1
                        du, dv = D*NB[:, k], D*NB[:, l]
                        duf(x) = SCeval(nodes, du, x)
                        dvf(x) = SCeval(nodes, dv, x)
                        L[k, l], err = quadgk(x->duf(x)*Visc(x)*dvf(x), -1.0, 1.0)
                end
        end
        S2 = -inv(Mref)*L
        w, v = eigen(-S2)
        locdt = 0.1/w[end]      # control the timestep using the eigenvalues
        disc(u, p, t) = S2*u
        p = 0.0
        prob = ODEProblem(disc, u0ar, tmax, p)
        sol = solve(prob, SSPRK33(), dt=locdt)
end


# Function calculates a set of dissipation matrices
function DispMats(xar, war)
        # calculation of the heat equation solution
        p = length(xar) - 1
	u0ar = Float64.(I[1:p+1, 1:p+1])
        sol = SolHESGal(u0ar, 10.0)
	# determining the minimal smoothing time for a positive solution
        # using bisection
        Di = zeros(p+1, p+1)
        tup = 10.0
        tlow = 0.0
        for it=1:16
                tmid = (tup + tlow) / 2
                Di[:] = sol(tmid)
                if sum(Di .< 0) > 0
			tlow = tmid
                else
                        tup = tmid
                end
        end
        println("Found positive heat equation solution after t= ", tup)
        # Collecting the matrix elements
        Di[:, :] = sol(tup)
        # check the magnitudes
        for k=1:p+1
                hges = war[k]
                hsol = dot(Di[:, k], war)
                println("Soll: ", hges, " Ist: ", hsol)
                hsol = sum(Di[k, :])
                println("Soll: ", 1.0, " Ist: ", hsol)

        end

        # Calcualting the generator
        G = Di - I

        println("opnorm ", opnorm(G))
        return G/max(tup, opnorm(G))
end

