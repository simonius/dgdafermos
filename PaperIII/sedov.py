# This file produces the sedov test pictures in the manuscript.

# Inclusion of the import module
from solimport import *
from schlieren import *

# Plots the grid in Pset and Tset, with collocation points pts and calues vals
# into files in figures/ starting with $Name.
def SedovPlot(name, Pset, Tset, pts, vals, p, ptpts, u):
    # Calculates the schlieren map
    smap = schlieren.ptschl(ptpts, u, 1)

    # Plot of the grid
    PlotGrid(Pset, Tset)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"grid.png", dpi=600)
    plt.clf()

    # A density map
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.colorbar()
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")

    plt.savefig("figures/"+name+"rho.png", dpi=600)
    plt.clf()

    # Schlieren images 
    plt.tricontourf(pts[:, 0], pts[:, 1], np.sqrt(smap), levels=1000, cmap='Greys')
    plt.title("Density, Schlieren image")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"smap.png", dpi=600)
    plt.clf()

    # pressure plots
    plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
    plt.colorbar()
    plt.title("Pressure")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"p.png", dpi=600)
    plt.clf()


Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/sedovp3gridexport.h5", "output/sedovp3solution.h5")
SedovPlot("sedovp3", Pset, Tset, pts, vals, p, ptpts, u)

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/sedovp1gridexport.h5", "output/sedovp1solution.h5")
SedovPlot("sedovp1", Pset, Tset, pts, vals, p, ptpts, u)

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/sedovp2gridexport.h5", "output/sedovp2solution.h5")
SedovPlot("sedovp2", Pset, Tset, pts, vals, p, ptpts, u)

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/sedovp5gridexport.h5", "output/sedovp5solution.h5")
SedovPlot("sedovp5", Pset, Tset, pts, vals, p, ptpts, u)

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/sedovp7gridexport.h5", "output/sedovp7solution.h5")
SedovPlot("sedovp7", Pset, Tset, pts, vals, p, ptpts, u)
