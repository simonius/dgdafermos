# File contains the Flux funcitons and approximate Riemann solvers for the 1D Euler equations
# u = (rho, rho u, E)

# the adiabatic exponent
const gamma = 1.4
# the number of conserved variables
const ncons = 3

# calculates the pressure from the conserved variables
function p(u)
	p = (gamma-1)*(u[3]-0.5*u[2]^2/u[1])
	return max(p, 1.0E-30)
end

# calculates the speed of sound
function a(u, p)
	a = sqrt(abs(gamma*p / u[1]))
        return a
end

# short form.
a(u) = a(u, p(u))

# The numerical speed of sound. Used instead of the exact sound speed in computations of al and ar
# to reduce the sonic point glitch, see also Tang: On the sonic point glitch
function anum(u, p)
        a = sqrt(abs(3*p / u[1]))
        return a
end

# the continuous Euler flux function
function fEuler(u, p)
        fEuler = zeros(3)
        fEuler[1] = u[2]
        fEuler[2] = u[2]^2/u[1] + p
        fEuler[3] = u[2]/u[1]*(u[3]+p)
        return fEuler
end

# The derivative of the Euler flux function with respect to $u$
function dfEuler(u)
	dfEuler = zeros(3, 3)
	dfEuler[1, :] = [0, 1, 0]
	dfEuler[2, :] = [-0.5*(gamma - 3)*(u[2]/u[1])^2, (3-gamma)*(u[2]/u[1]), gamma-1]
	dfEuler[3, 1] = - gamma*u[2]*u[3]*u[1]^(-2) + (gamma-1)*(u[2]/u[1])^3
	dfEuler[3, 2] = gamma*u[3]/u[1] - 1.5*(gamma-1)*(u[2]/u[1])^2 
	dfEuler[3, 3] = gamma*u[2]/u[1]
	return dfEuler
end

# short form for fEuler(u, p)
function fEuler(u)
        return fEuler(u, p(u))
end


fopnorm(ul, ur) = max(abs(ul[2]/ul[1]), abs(ur[2]/ur[1])) + max(anum(ul, p(ul)), anum(ur, p(ur)))
fsb(ul, ur) = (min(ul[2] / ul[1], ur[2] / ur[1]) - max(anum(ul, p(ul)), anum(ur, p(ur))), max(ul[2] / ul[1], ur[2] / ur[1]) + max(anum(ul, p(ul)), anum(ur, p(ur))))


# Lax-Friedrichs flux for the Euler equations
function hEulerLF(ul, ur, mu)
        hEulerLF = zeros(3)
        pl = p(ul)
        pr = p(ur)
        hEulerLF = 0.5*(fEuler(ul, pl) + fEuler(ur, pr) - (ur-ul)*mu)
        return hEulerLF
end

# Entropy variables (derivative of the entropy functional) for the Euler equations
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


# Local Lax-Friedrichs flux for the Euler equations
@inline function hEulerLLF(ul, ur)
        pl = p(ul)
        pr = p(ur)
        amax = max(anum(ul, pl), anum(ur, pr))
	lambound = fopnorm(ul, ur)
	return hEulerLF(ul, ur, lambound)
end

# Local Lax-Friedrichs entropy flux for the Euler equations
@inline function HEulerLLF(ul, ur)
        pl = p(ul)
        pr = p(ur)
        amax = max(anum(ul, pl), anum(ur, pr))
	lambound = fopnorm(ul, ur)
	return HEulerLF(ul, ur, lambound)
end


# The physical Euler Entropy see Harten (1983) On the symmetric form of systems of conservation laws
function PUEuler(u)
        S = log(p(u)*abs(u[1])^(-gamma))
        return -u[1]*S
end

# The entropy flux for the physical entropy and the Euler equations
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

# The entropy flux for the Lax Friedrichs flux
function HEulerLF(ul, ur, mu)
        return 0.5*(PFEuler(ul) + PFEuler(ur) - (PUEuler(ur) - PUEuler(ul))*mu)
end

# The HLL numerical flux for the Euler equations of gas dynamics
function hEulerHLL(ul, ur)
        hEulerHLL = zeros(3)
        pl = p(ul)
        pr = p(ur)
        amax = max(anum(ul, pl), anum(ur, pr))
        vl = ul[2]/ul[1]
        vr = ur[2]/ur[1]
        vmin = min(vl, vr)
        vmax = max(vl, vr)
        al = vmin-amax
        ar = vmax+amax

        if (0 < al)
                hEulerHLL = fEuler(ul,pl)
        elseif (0 < ar)
                hEulerHLL = (al*ar*(ur-ul) + ar*fEuler(ul, pl) - al*fEuler(ur, pr)) / (ar-al)
        else
                hEulerHLL = fEuler(ur, pr)
        end
        return hEulerHLL
end

# The HLL numerical entropy flux for the Euler equations of gas dynamics
function HEulerHLL(ul, ur)
        HEulerHLL = zeros(3)
        pl = p(ul)
        pr = p(ur)
        amax = max(anum(ul, pl), anum(ur, pr))
        vl = ul[2]/ul[1]
        vr = ur[2]/ur[1]
        vmin = min(vl, vr)
        vmax = max(vl, vr)
        al = vmin-amax
        ar = vmax+amax
        um = um = (ar*ur - al*ul)/(ar-al) + (fEuler(ul) - fEuler(ur))/(ar-al)

        if (0 < al)
                flb = al*PUEuler(ul) + (ar - al) * PUEuler(um) - ar*PUEuler(ur) + PFEuler(ur)
                fub = PFEuler(ul)
                HEulerHLL = fub
        elseif (0 < ar)
                flb = ar*PUEuler(um) - ar*PUEuler(ur) + PFEuler(ur)
                fub = al*PUEuler(um) - al*PUEuler(ul) + PFEuler(ul)
         #       HEulerHLL = -al/(ar - al)*flb + ar/(ar - al)*fub
                HEulerHLL = 0.5*(flb + fub)
        else
                flb = PFEuler(ur)
                fub = PFEuler(ul) - (ar - al)*PUEuler(um)+ ar*PUEuler(ur) - al*PUEuler(ul)
                HEulerHLL = flb
        end
        return HEulerHLL
end

# The Dissipation Prediction of a local Lax-Friedrichs style entropy inequality predictor,
# i.e. the HLL predictor from the publication with the trivial speed estimate -a_l = c = a_r
function LLFDispOp(ul, ur)
        pl, pr = p(ul), p(ur)
        amax = max(anum(ul, pl), anum(ur, pr))
        cbound = fopnorm(ul, ur)
        um = 0.5*(ul + ur) + (fEuler(ul) - fEuler(ur))/(2*cbound)
        return cbound*(2*Ufunc(um) - Ufunc(ul) - Ufunc(ur)) - F(ul) + F(ur)
end

# The Dissipation prediction of a HLL style entropy inequality predictor.
function HLLDispOp(ul, ur)
        pl, pr = p(ul), p(ur)
        amax = max(anum(ul, pl), anum(ur, pr))
        ar = max(ul[2]/ul[1], ur[2]/ur[1]) + amax
        al = min(ul[2]/ul[1], ur[2]/ur[1]) - amax
        um = (ar*ur - al*ul)/(ar-al) + (fEuler(ul) - fEuler(ur))/(ar-al)
        lUd = (ar - al)*Ufunc(um) + al*Ufunc(ul) - ar*Ufunc(ur) - (F(ul) - F(ur))

        # We do not know in which cell the dissipation takes place
        return lUd
end

# Returns the conserved variables given the physical variables density, speed and pressure
function cons_var(rho, v, p)
    E = p / (gamma-1) + 0.5*v^2*rho
    rhov = rho*v
    return [rho, rhov, E]
end

# Initial conditions
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


function toro1(x)     
        if x < 3.0
                u0 = cons_var(1.0, 0.75, 1.0)
        else
                u0 = cons_var(0.125, 0.0, 0.10)
        end

        return u0
end

function toro3(x)
        if x < 5.0
                u0 = cons_var(1.0, 0.0, 1000.0)
        else
                u0 = cons_var(1.0, 0.0, 0.01)
        end

        return u0
end

function toro4(x)
        if x < 4.0
                u0 = cons_var(5.99924, 19.5975, 460.894)
        else
                u0 = cons_var(5.99242, -6.19633, 46.0950)
        end

        return u0
end

function toro5(x)
        if x < 8.0
                u0 = cons_var(1.0, -19.59745, 1000.0)
        else
                u0 = cons_var(1.0, -19.59745, 0.01)
        end

        return u0
end


function smootheuler(x)
        eps = 0.2
        
        eps = exp(-(x-5.0)^2)
        u0 = cons_var(3.85 + eps*sin(2*x), 2.0, 10.3)
        
        return u0
end


