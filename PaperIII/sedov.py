# This file produces the sedov test pictures in the manuscript.

# Inclusion of the import module
from solimport import *

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/sedovgridexport.h5", "output/sedovsolution.h5")

# Plots the grid in Pset and Tset, with collocation points pts and calues vals
# into files in figures/ starting with $Name.
def SedovPlot(name, Pset, Tset, pts, vals, p, Maz):
    PlotGrid(Pset, Tset)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"grid.png", dpi=600)
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.colorbar()
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"rho.png", dpi=600)
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
    plt.colorbar()
    plt.title("Pressure")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"p.png", dpi=600)
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], Maz, levels=1000)
    plt.colorbar()
    plt.title("Mach Number")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"M.png", dpi=600)
    plt.clf()

SedovPlot("sedov", Pset, Tset, pts, vals, p, Maz)

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/sedovp1gridexport.h5", "output/sedovp1solution.h5")
SedovPlot("sedovp1", Pset, Tset, pts, vals, p, Maz)
