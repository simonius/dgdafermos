using DifferentialEquations

include("euler.jl")


function EulerLLF(du, u, p, t)
	dx = p[1]
	dxi = 1/dx
	
	N::Int64 = size(u)[1]
	
	Threads.@threads for i = 3:N-2
                du[i, :] = (hEulerLLF(u[i-1, :], u[i, :]) - hEulerLLF(u[i, :], u[i+1, :]))*dxi
        end
        du[1:2, :] .= 0
       	du[N-1:N, :] .= 0
end

function EulerLF(u0f, tmax, N, cbound)
	fa = zeros(N, 3)
	dx = (xr - xl) / N
	dt = dx/cbound
	lamh = 0.5*dt/dx
	K = ceil(Int, tmax/dt)
	u = zeros(K, N, 3)
	Threads.@threads for k=1:N
                u[1, k, :] = u0f(xl + (k-0.5)*dx)
        end

	for k=1:K-1
		Threads.@threads for n=1:N
			fa[n, :] = fEuler(u[k, n, :], p(u[k, n, :]))
		end
		u[k+1, 2:end-1,:] = 0.5*(u[k, 1:end-2, :] + u[k, 3:end, :]) + lamh*(fa[1:end-2, :] - fa[3:end, :])
	
		u[k+1, 1, :] = u[k, 1, :]
		u[k+1, end, :] = u[k, end, :]
	end
	return u
end

# Same as before, but only saves final value
function EulerLFlu(u0f, tmax, N, cbound)
        fa = zeros(N, 3)
        dx = (xr - xl) / N
        dt = dx/cbound
        lamh = 0.5*dt/dx
	K = ceil(Int, tmax/(2*dt))
        u = zeros(N, 3)
	u2 = zeros(N, 3)
        Threads.@threads for k=1:N
                u[k, :] = u0f(xl + (k-0.5)*dx)
        end

        for k=1:K
		# u->u2
                Threads.@threads for n=1:N
                        fa[n, :] = fEuler(u[n, :], p(u[n, :]))
                end
                u2[2:end-1,:] = 0.5*(u[1:end-2, :] + u[3:end, :]) + lamh*(fa[1:end-2, :] - fa[3:end, :])

                u2[1, :] = u[1, :]
                u2[end, :] = u[end, :]

		# u2->u1
		Threads.@threads for n=1:N
                        fa[n, :] = fEuler(u2[n, :], p(u2[n, :]))
                end
                u[2:end-1,:] = 0.5*(u2[1:end-2, :] + u2[3:end, :]) + lamh*(fa[1:end-2, :] - fa[3:end, :])

                u[ 1, :] = u2[1, :]
                u[end, :] = u2[end, :]
        end
        return u
end


function CalcRefNPFV(solver;u0f = shuosh6, tmax=1.8, N=100, CFL=1.0, cbound = 3.0)
	dx = (xr - xl)/N
	u0ar = zeros(N, 3)
	for k=1:N
		u0ar[k, :] = u0f(xl + (k-0.5)*dx)
	end
        p = zeros(2)
        p[1] = dx
        locdt = CFL*p[1]/cbound
	p[2] = locdt
        prob = ODEProblem(solver, u0ar, tmax, p)
	# We are using SSPRK33 as SSPRK22 has bad eigenvalues (Private Comm with Kurganov)
        sol = solve(prob, SSPRK33(), dt=locdt,  saveat=0.01)
        return sol
end

