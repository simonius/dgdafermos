# This file implements the needed basis changes and inner product
# Grammians for a classical DG method with the exact Mass matrix.
using LinearAlgebra
using QuadGK

# recursively gives the Lagrange polynomials
# at point x, counting from 0
function Pleg(n, x)
	if n == 0
		return 1.0
	elseif n == 1
		return x
	else
		return ((2*n-1)*x*Pleg(n-1, x) - (n-1)*Pleg(n-2, x))/n
	end
end

# recursively gives the derivative of
# Lagrange polynomial $n$, counting from zero,
# at position x
function dPleg(n, x)
	if n == 0
		return 0.0
	elseif n == 1
		return 1.0
	else
		return ((2*n-1)*(Pleg(n-1, x) + x*dPleg(n-1, x)) - (n-1)*dPleg(n-2, x))/n
	end
end

# Matrix entries of the modal mass matrix
# in the Lagrange Basis, obv. diagonal
# counting from 0
function Me(k, l)
	if k == l
		return 2/(2*k+1)
	else
		return 0.0
	end
end

# Matrix entries of the stiffness matrix
# for the first derivative
# counting from 0
function Se(k, l)
	if l == 0
		return 0.0
	elseif l == 1
		return Me(k, 0)
	else
		return (2*l-1)*Me(k, l-1) + Se(k, l-2)
	end
end

# Helper function, filling a matrix with 
# values from one of the above functions
function Mat(ME, N)
	M = zeros(N, N)
	for k=1:N
		for l=1:N
			M[k, l] = ME(k-1, l-1)
		end
	end
	return M
end

# Fills a matrix with columns of at x 
# evaluated Lagrange polynomials, 
# allowing to change from modal
# to nodal basis.
function PaMat(x)
	N = length(x)
	M = zeros(N, N)
	for k=1:N
		for l=1:N
			M[k, l] = Pleg(l-1, x[k])
		end
	end
	return M
end


# calculates the Lagrange Polynomial k to x_k at x
function Plag(xar, k, x)
        K = length(xar)
        v = 1.0
        for i=1:k-1
                v = v*(x-xar[i])/(xar[k] - xar[i])
        end
        for i=k+1:K
                v = v*(x-xar[i])/(xar[k] - xar[i])
        end
        return v
end

# calcualtes the value of $u$ at $x$
# given a vector of nodal coefficients $uT$
function SCeval(xar, uT, x)
        v = 0.0
        K = length(xar)
        for k=1:K
                v = v + Plag(xar, k, x)*uT[k]
        end
        return v
end

# calculates the value of $u$ at $x$
# given a vector of Legendre modal coefficients c
function CellEval(c, x)
	N = length(c)
	s = 0.0
	for n=1:N
		s = s + c[n]*Pleg(n-1, x)
	end
	return s
end
