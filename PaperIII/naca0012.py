# This file is used to produce all graphs for the NACA 0012 testcase
# in the manuscript.

from solimport import *
from matplotlib import patches
from schlieren import schlieren

# import the airfoil for solution overlay
af = np.genfromtxt("grids/Naca0012", skip_header=1)
af[:, 0] = af[:, 0] - 0.5

def NacaPlot(name, af, Pset, Tset, pts, vals, p, ptpts, u):
     # Calculates the schlieren map
    smap = schlieren.ptschl(ptpts, u, 1)

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
    airfoil = patches.Polygon(af, facecolor='white', zorder = 10)
    ax.add_patch(airfoil)

    plt.savefig("figures/"+name+"rho.png", dpi=600,bbox_inches='tight')
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.savefig("figures/"+name+"rhoclose.png", dpi=600,bbox_inches='tight')
    plt.clf()

    # Schlieren images 
    plt.tricontourf(pts[:, 0], pts[:, 1], np.sqrt(smap), levels=1000, cmap='Greys')
    plt.title("Density, Schlieren image")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    ax = plt.gca()
    airfoil = patches.Polygon(af, facecolor='white', zorder = 10)
    ax.add_patch(airfoil)
    plt.savefig("figures/"+name+"smap.png", dpi=600,bbox_inches='tight')
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.savefig("figures/"+name+"smapclose.png", dpi=600,bbox_inches='tight')
    plt.clf()

    # Density with grid overlay
    PlotGrid(Pset, Tset)
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.colorbar()
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    ax = plt.gca()
    airfoil = patches.Polygon(af, facecolor='white', zorder = 10)
    ax.add_patch(airfoil)
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.savefig("figures/"+name+"rhogrid.png", dpi=600,bbox_inches='tight')
    plt.clf()

    # Pressure
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
 

# Solutions on a fine grid up to $p=3$
Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/Naca0012p1gridexport.h5", "output/Naca0012p1solution.h5")
NacaPlot("Nacap1", af, Pset, Tset, pts, vals, p, ptpt, u)

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/Naca0012p2gridexport.h5", "output/Naca0012p2solution.h5")
NacaPlot("Nacap2", af, Pset, Tset, pts, vals, p, ptpt, u)

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/Naca0012p3gridexport.h5", "output/Naca0012p3solution.h5")
NacaPlot("Nacap3", af, Pset, Tset, pts, vals, p, ptpt, u)

# Solutions on a coarse grid from $p=3$
Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/NacaLRp3gridexport.h5", "output/NacaLRp3solution.h5")
NacaPlot("NacaLRp3", af, Pset, Tset, pts, vals, p, ptpt, u)

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/NacaLRp5gridexport.h5", "output/NacaLRp5solution.h5")
NacaPlot("NacaLRp5", af, Pset, Tset, pts, vals, p, ptpt, u)

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/NacaLRp7gridexport.h5", "output/NacaLRp7solution.h5")
NacaPlot("NacaLRp7", af, Pset, Tset, pts, vals, p, ptpt, u)


