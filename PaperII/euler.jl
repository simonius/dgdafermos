# File contains the Flux funcitons and Riemann solvers for the 1D Euler equations
# u = (rho, rho u, E)

gamma = 1.4
ncons = 3

function p(u)
	p = (gamma-1)*(u[3]-0.5*u[2]^2/u[1])
	return max(p, 1.0E-30)
end

function a(u, p)
	a = sqrt(abs(gamma*p / u[1]))
        return a
end

function fEuler(u, p)
        fEuler = zeros(3)
        fEuler[1] = u[2]
        fEuler[2] = u[2]^2/u[1] + p
        fEuler[3] = u[2]/u[1]*(u[3]+p)
        return fEuler
end

function dfEuler(u)
	dfEuler = zeros(3, 3)
	dfEuler[1, :] = [0, 1, 0]
	dfEuler[2, :] = [-0.5*(gamma - 3)*(u[2]/u[1])^2, (3-gamma)*(u[2]/u[1]), gamma-1]
	dfEuler[3, 1] = - gamma*u[2]*u[3]*u[1]^(-2) + (gamma-1)*(u[2]/u[1])^3
	dfEuler[3, 2] = gamma*u[3]/u[1] - 1.5*(gamma-1)*(u[2]/u[1])^2 
	dfEuler[3, 3] = gamma*u[2]/u[1]
	return dfEuler
end

function fEuler(u)
        return fEuler(u, p(u))
end

fopnorm(ul, ur) = max(abs(ul[2]/ul[1] - a(ul, p(ul))), abs(ur[2]/ur[1] + a(ur, p(ur))))


function hEulerLF(ul, ur, mu)
        hEulerLF = zeros(3)
        pl = p(ul)
        pr = p(ur)
        hEulerLF = 0.5*(fEuler(ul, pl) + fEuler(ur, pr) - (ur-ul)*mu)
        return hEulerLF
end


function eulerEvar(cu)
        rho = cu[1]
        u = cu[2]/cu[1]
        E = cu[3]
        p = (gamma-1)*(E - 0.5*rho*u^2)
        s = log(p) - gamma*log(rho)
        v1 = (gamma-s)/(gamma-1) - (rho*u^2) / (2*p)
        v2 = rho*u/p
        v3 = - rho / p
        return [v1, v2, v3]
end



@inline function hEulerLLF(ul, ur)
        pl = p(ul)
        pr = p(ur)
        amax = max(a(ul, pl), a(ur, pr))
	lambound=  1.0*(max(abs(ul[2]/ul[1]), abs(ur[2]/ur[1])) + amax)
        return hEulerLF(ul, ur, lambound)
end

@inline function HEulerLLF(ul, ur)
        pl = p(ul)
        pr = p(ur)
        amax = max(a(ul, pl), a(ur, pr))
	lambound=  1.0*(max(abs(ul[2]/ul[1]), abs(ur[2]/ur[1])) + amax)
        return HEulerLF(ul, ur, lambound)
end


# The physical Euler Entropy see Harten (1983) On the symmetric form of systems of conservation laws
function PUEuler(u)
        S = log(p(u)*abs(u[1])^(-gamma))
        return -u[1]*S
end

function PFEuler(u)
        S = log(p(u)*abs(u[1])^(-gamma))
        return -u[2]*S
end

# Entropy variables, Harten (1983b)/ Tadmor 2003 Acta Numerica
function dUEuler(u)
	S = log(abs(p(u)*abs(u[1])^(-gamma)))
	h(S) = S
	h´(S) = 1.0
	v = (1-gamma)*h´(S)/p(u)*[u[3] + p(u)/(gamma-1)*(h(S)/h´(S) - gamma - 1),
				  -u[2],
				  u[1]]
	return v
end


function HEulerLF(ul, ur, mu)
        return 0.5*(PFEuler(ul) + PFEuler(ur) - (PUEuler(ur) - PUEuler(ul))*mu)
end


function cons_var(rho, v, p)
    E = p / (gamma-1) + 0.5*v^2*rho
    rhov = rho*v
    return [rho, rhov, E]
end


function shuosh6(x)
        eps = 0.2
        
        if x < 1.0
                u0 = cons_var(3.857143, 2.629369, 10.33333333333)
        else
                u0 = cons_var(1.0 + eps*sin(5*x), 0.0, 1.0)
        end
        
        return u0
end

function shuosh6a(x)
        
        if x < 5.0
                u0 = cons_var(1.0, 0.0, 1.0)
        else
                u0 = cons_var(0.125, 0.0, 0.10)
        end
       
        return u0
end

function shuosh6b(x)
        if x < 5.0
                u0 = cons_var(0.445, 0.698, 3.528)
        else
                u0 = cons_var(0.5, 0.0, 0.571)
        end
        return u0
end


function shuosh7(x)
        ul = cons_var(1.0, 0.0, 1E3)
        um = cons_var(1.0, 0.0, 1E-2)
        ur = cons_var(1.0, 0.0, 1E2)
        if x < 1.0
                u0 = ul
        elseif x < 9.0
                u0 = um
        else
                u0 = ur
        end
        return u0
end


function toro123(x)
        if x < 5.0
                u0 = cons_var(1.0, -2.0, 0.4)
        else
                u0 = cons_var(1.0, 2.0, 0.4)
        end
        return u0
end



function smootheuler(x)
        eps = 0.2
        
        eps = exp(-(x-5.0)^2)
        u0 = cons_var(3.85 + eps*sin(2*x), 2.0, 10.3)
        
        return u0
end


