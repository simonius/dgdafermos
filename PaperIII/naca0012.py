# This file is used to produce all graphs for the NACA 0012 testcase
# in the manuscript.

from solimport import *
from matplotlib import patches

# import the airfoil for solution overlay
af = np.genfromtxt("grids/Naca0012", skip_header=1)
af[:, 0] = af[:, 0] - 0.5
Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/Naca0012gridexport.h5", "output/Naca0012solution.h5")

def NacaPlot(name, af, Pset, Tset, pts, vals, p, Maz):
    PlotGrid(Pset, Tset)
    plt.axis("equal")
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("figures/"+name+"grid.png", dpi=600,bbox_inches='tight')
    plt.xlim(-0.6, 0.6)
    plt.ylim(-0.6, 0.6)
    plt.savefig("figures/"+name+"gridclose.png", dpi=600,bbox_inches='tight')
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.colorbar()
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    ax = plt.gca()
    airfoil = patches.Polygon(af, facecolor='white')
    ax.add_patch(airfoil)
    plt.savefig("figures/"+name+"rho.png", dpi=600,bbox_inches='tight')
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.savefig("figures/"+name+"rhoclose.png", dpi=600,bbox_inches='tight')
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
    plt.colorbar()
    plt.title("Pressure")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    ax = plt.gca()
    airfoil = patches.Polygon(af, facecolor='white')
    ax.add_patch(airfoil)
    plt.savefig("figures/"+name+"p.png", dpi=600,bbox_inches='tight')
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.savefig("figures/"+name+"pclose.png", dpi=600,bbox_inches='tight')
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], Maz, levels=1000)
    plt.colorbar()
    plt.title("Mach Number")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    ax = plt.gca()
    airfoil = patches.Polygon(af, facecolor='white')
    ax.add_patch(airfoil)
    plt.savefig("figures/"+name+"M.png", dpi=600,bbox_inches='tight')
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.savefig("figures/"+name+"Mclose.png", dpi=600,bbox_inches='tight')
    plt.clf()



NacaPlot("Naca", af, Pset, Tset, pts, vals, p, Maz)

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/Naca0012M2gridexport.h5", "output/Naca0012M2solution.h5")
NacaPlot("NacaM2", af, Pset, Tset, pts, vals, p, Maz)


Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/Naca0012p1gridexport.h5", "output/Naca0012p1solution.h5")
NacaPlot("Nacap1", af, Pset, Tset, pts, vals, p, Maz)

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/Naca0012p1M2gridexport.h5", "output/Naca0012p1M2solution.h5")
NacaPlot("Nacap1M2", af, Pset, Tset, pts, vals, p, Maz)


Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/NacaLRgridexport.h5", "output/NacaLRsolution.h5")
NacaPlot("NacaLR", af, Pset, Tset, pts, vals, p, Maz)

