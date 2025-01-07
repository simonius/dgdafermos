# This file produces the plots of the forward facing step problem
# from the saved HDF5 files in output

from solimport import *
from matplotlib import patches
from schlieren import schlieren

def ffsplot(Pset, Tset, pts, vals, p, ptpts, u, name):
    # Calculates the schlieren map
    smap = schlieren.ptschl(ptpts, u, 1)


    # A plot of the grid
    PlotGrid(Pset, Tset)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("figures/"+name+"grid.png", dpi=600, bbox_inches='tight')
    plt.clf()

    # A plot of the density
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    rect = patches.Rectangle((0.6,0.0),2.4,0.2, facecolor='white')
    ax.add_patch(rect)
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)

    plt.savefig("figures/"+name+"rho.png", dpi=600, bbox_inches='tight')
    plt.clf()

    # Schlieren images
    plt.tricontourf(pts[:, 0], pts[:, 1], np.sqrt(smap), levels=1000, cmap='Greys')
    plt.title("Density, Schlieren image")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    rect = patches.Rectangle((0.6,0.0),2.4,0.2, facecolor='white')
    ax.add_patch(rect)
    plt.savefig("figures/"+name+"smap.png", dpi=600)
    plt.clf()


    # A plot of the pressure
    plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
    plt.title("Pressure")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    rect = patches.Rectangle((0.6,0.0),2.4,0.2, facecolor='white')
    ax.add_patch(rect)
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)
    plt.savefig("figures/"+name+"p.png", dpi=600, bbox_inches='tight')
    plt.clf()

# Tests for end time t = 0.5
Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs05p1gridexport.h5", "output/ffs05p1solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs05p1")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs05p2gridexport.h5", "output/ffs05p2solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs05p2")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs05p3gridexport.h5", "output/ffs05p3solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs05p3")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs05p5gridexport.h5", "output/ffs05p5solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs05p5")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs05p7gridexport.h5", "output/ffs05p7solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs05p7")

# Tests for end time t = 3.0
Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs30p1gridexport.h5", "output/ffs30p1solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs30p1")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs30p2gridexport.h5", "output/ffs30p2solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs30p2")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs30p3gridexport.h5", "output/ffs30p3solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs30p3")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs30p5gridexport.h5", "output/ffs30p5solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs30p5")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/ffs30p7gridexport.h5", "output/ffs30p7solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "ffs30p7")

