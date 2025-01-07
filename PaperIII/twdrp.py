# This file produces the plots of the forward facing step problem
# from the saved HDF5 files in output

from solimport import *
from matplotlib import patches
from schlieren import schlieren


def ffsplot(Pset, Tset, pts, vals, p, ptpts, u, name):
    # Calculates the schlieren map
    smap = schlieren.ptschl(ptpts, u, 1)

    # Grid plot
    PlotGrid(Pset, Tset)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("figures/"+name+"grid.png", dpi=600, bbox_inches='tight')
    plt.clf()

    # Density plot
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
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
    plt.savefig("figures/"+name+"smap.png", dpi=600)
    plt.clf()


    # Pressure plot
    plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
    plt.title("Pressure")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)
    plt.savefig("figures/"+name+"p.png", dpi=600, bbox_inches='tight')
    plt.clf()



Pset, Tset, pts, vals, p, ptpt, u= SolImport("output/twdrp1gridexport.h5", "output/twdrp1solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpt, u, "twdrp1")

Pset, Tset, pts, vals, p, ptpt, u= SolImport("output/twdrp2gridexport.h5", "output/twdrp2solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpt, u, "twdrp2")

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/twdrp3gridexport.h5", "output/twdrp3solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpt, u, "twdrp3")

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/twdrp5gridexport.h5", "output/twdrp5solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpt, u, "twdrp5")

Pset, Tset, pts, vals, p, ptpt, u = SolImport("output/twdrp7gridexport.h5", "output/twdrp7solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpt, u,  "twdrp7")

