# File that installs all needed packages
import Pkg
Pkg.add("FastGaussQuadrature")
Pkg.add("QuadGK")

Pkg.add("DifferentialEquations")
Pkg.add("Plots")
Pkg.add("Makie")
Pkg.add("CairoMakie")
Pkg.add("GLMakie")
Pkg.add("LaTeXStrings")
Pkg.update()

