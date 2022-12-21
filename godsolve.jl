# This file implements the reference solution solver used for shock capturing calculations
# the Godunov flux for Burgers equation is defined in the file for the DG solvers.

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
        u0 = u0func.(collect(grid))
        tspan = (0, tmax)
        prob = ODEProblem(BurgersGod!, u0, tspan, [dx, order])
        sol = solve(prob, Euler(), dt=dt)
        return sol
end

