# An exact Riemann Solver, based on Toros book.

# inclusion of the continuous Euler flux and formuals for sound speed and pressure
include("euler.jl")

# the functions fl and fr from toros book, depending on the 
# left and right state and the pressure in the star region
function fl(pstar, ul)
	if pstar > p(ul)
		Al = 2 / ((gamma +1) * ul[1])
		Bl = (gamma - 1) / (gamma + 1)*p(ul)
		return (pstar - p(ul))*sqrt(Al / (pstar + Bl))
	else
		return 2.0*a(ul)/(gamma-1.0)*((pstar / p(ul))^((gamma-1.0)/(2.0*gamma)) - 1.0)
	end
end

function fr(pstar, ur)
        if pstar > p(ur)
		Ar = 2 / ((gamma +1) * ur[1])
                Br = (gamma - 1) / (gamma + 1)*p(ur)
		return (pstar - p(ur))*sqrt(Ar / (pstar + Br))
        else    
		return 2.0*a(ur)/(gamma-1.0)*((pstar / p(ur))^((gamma-1.0)/(2.0*gamma)) - 1.0)
        end     
end

# the combined function whose zero yields pstart
function pdiskr(pstar, ul, ur)
	vl, vr = ul[2]/ul[1], ur[2]/ur[1]
	return fl(pstar, ul) + fr(pstar, ur) + vr - vl
end

# this subroutine determines pstar from ul and ur
# Niter is the number of bisections used to find the zero of pdiskr
# 64 iterations, i.e. 64 bits of precision should be enough for double precision.
function pstar(ul, ur, Niter=64)
	plb = 0.0
	# first step: search for an upper bound of p
	pub = 2
	for k=1:10
		if pdiskr(pub, ul, ur) < 0.0
			plb = pub
			pub = 2*pub
		end
	end
	# now bisect onto the pressure in the star region
	for k=1:Niter
		if pdiskr(0.5*(plb + pub), ul, ur) < 0.0
			plb = 0.5*(plb + pub)
		else
			pub = 0.5*(plb + pub)
		end
	end
	return 0.5*(plb+pub)
end

# calculates the speed in the star region
function vstar(ps, ul, ur)
	vl, vr = ul[2]/ul[1], ur[2]/ur[1]
	return 0.5*(vl + vr + fr(ps, ur) - fl(ps, ul))
end

# returns the density left of the contact if the left wave is a shock
# and the speed of the left going shock wave
function ShockL(ps, vs, ul)
	pl = p(ul)
	rhosl = ul[1]*(ps / pl + (gamma-1)/(gamma+1)) / ((gamma-1) / (gamma+1) * ps/pl + 1)
        Al = 2 / ((gamma +1) * ul[1])
        Bl = (gamma - 1) / (gamma + 1)*p(ul)
        Ql = sqrt((ps + Bl) / Al)
        Sl = ul[2]/ul[1] - Ql / ul[1]
        return rhosl, Sl
end

# returns the density left of the contact and 
# the head and tail speed if the left wave is a rarefaction
function RarefacL(ps, vs, ul)
	rhosl = ul[1]*(ps/p(ul))^(1/gamma)
        asl = a(ul)*(ps/p(ul))^((gamma-1)/(2*gamma))
        Shl = ul[2]/ul[1] - a(ul)
        Stl = vs - asl
        rhosl, Shl, Stl
end

# value of the conserved quantities in a left rarefaction fan
function LFan(ul, xbt)
        vl = ul[2] / ul[1]
        rhofan = ul[1]*(2/(gamma+1) + (gamma-1)/((gamma+1)*a(ul))*(vl - xbt))^(2/(gamma-1))
        vfan = 2 / (gamma + 1) * (a(ul) + (gamma - 1) / 2 * vl + xbt)
        pfan = p(ul)*(2 / (gamma + 1) + (gamma-1) / ((gamma + 1) * a(ul))*(vl - xbt) )^(2*gamma/(gamma-1))
        return cons_var(rhofan, vfan, pfan)
end

# returns the density right of the contact if the right wave is a shock
# and the speed of the right going shock wave
function ShockR(ps, vs, ur)
        pr = p(ur)
        rhosr = ur[1]*(ps / pr + (gamma-1)/(gamma+1)) / ((gamma-1) / (gamma+1) * ps/pr + 1)  
        #Qr = rhosr*vs
	Ar = 2 / ((gamma +1) * ur[1])
        Br = (gamma - 1) / (gamma + 1)*p(ur)
	Qr = sqrt((ps + Br) / Ar)
        Sr = ur[2]/ur[1] + Qr / ur[1]
        return rhosr, Sr
end

# returns the density right of the contact and
# the head and tail speed if the right wave is a rarefaction
function RarefacR(ps, vs, ur)
	rhosr = ur[1]*(ps/p(ur))^(1/gamma)
        asr = a(ur)*(ps/p(ur))^((gamma-1)/(2*gamma))
        Shr = ur[2]/ur[1] + a(ur)
        Str = vs + asr
        return rhosr, Str, Shr
end

# value of the conserved quantities in a right rarefaction fan
function RFan(ur, xbt)
	vr = ur[2] / ur[1]
        rhofan = ur[1]*(2/(gamma+1) - (gamma-1)/((gamma+1)*a(ur))*(vr - xbt))^(2/(gamma-1))
        vfan = 2 / (gamma + 1) * (-a(ur) + (gamma - 1) / 2 * vr + xbt)
        pfan = p(ur)*(2 / (gamma + 1) - (gamma-1) / ((gamma + 1) * a(ur))*(vr - xbt) )^(2*gamma/(gamma-1))
        return cons_var(rhofan, vfan, pfan)
end

# solves the Riemann problem for the initial values ul, ur
function uRiemann(ul, ur, xbt; debug = false)
	# determine star region
	ps = pstar(ul, ur)
	vs = vstar(ps, ul, ur)
	if debug
		println("star pressure is: ", ps)
		println("star speed is: ", vs)
		if (ps < p(ul))
			println("left rarefaction")
		else
			println("left shock")
		end
		if (ps < p(ur))
			println("right rarefaction")
		else
			println("right shock")
		end
	end
	# Now sample the solution
	if xbt < vs # left of contact
		if ps < p(ul) 	# left rarefaction
			rhosl, Shl, Stl = RarefacL(ps, vs, ul)
			if xbt < Shl
				return ul
			elseif xbt < Stl
				return LFan(ul, xbt)
			else
				return cons_var(rhosl, vs, ps)
			end
		else 		# left shock
			rhosl, Sl = ShockL(ps, vs, ul)
			if xbt < Sl
				return ul
			else
				return cons_var(rhosl, vs, ps)
			end
		end
	else
		if ps < p(ur)	# right rarefaction
			rhosr, Str, Shr = RarefacR(ps, vs, ur)
			if xbt < Str
				return cons_var(rhosr, vr, ps)
			elseif xbt < Shr
				return RFan(ur, xbt)
			else
				return ur
			end
		else		# right shock
			rhosr, Sr = ShockR(ps, vs, ur)
			if xbt > Sr
				return ur
			else
				return cons_var(rhosr, vs, ps)
			end
		end
	end
end
