# This file implements the reference solution solver used for shock capturing calculations 
# of Burgers' equation.


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


# ODE system
function BurgersGod!(du, u, p, t)
        dx = p[1]
        k = p[2]
        N = length(u)
        for i = 1:length(u)
                du[i] = (hBurgersGod(u[mod1(i-1, N)], u[i]) - hBurgersGod(u[i], u[mod1(i+1, N)]))/dx
        end
end

# Driver to solver this problem with forward Euler time-integration
function GS(N, tmax, u0func, xmax, cfl)
        #dx = xmax / N
        grid = (0.0:2/N:2.0)[1:N]
        dx = grid[2]-grid[1]
        dt = cfl*dx/4.0
        println("dx = " * string(dx))
        println("dt = " * string(dt))
        u0 = zeros(N)
        for n=1:N
                u0[n] = u0func(grid[n])[1]
        end
        tspan = (0, tmax)
        prob = ODEProblem(BurgersGod!, u0, tspan, [dx, order])
	sol = solve(prob, Euler(), dt=dt, saveat=[0.0, tmax])
        return sol
end

