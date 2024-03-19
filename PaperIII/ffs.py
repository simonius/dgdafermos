# This file produces the plots of the forward facing step problem
# from the saved HDF5 files in output

from solimport import *
from matplotlib import patches


def ffsplot(Pset, Tset, pts, vals, p, Maz, cb, name):
    PlotGrid(Pset, Tset)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("figures/"+name+"grid.png", dpi=600, bbox_inches='tight')
    plt.clf()

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

    plt.tricontourf(pts[:, 0], pts[:, 1], Maz, levels=1000, vmax=5.0)
    plt.title("Mach Number")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    rect = patches.Rectangle((0.6,0.0),2.4,0.2, facecolor='white')
    ax.add_patch(rect)
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)
    plt.savefig("figures/"+name+"M.png", dpi=600, bbox_inches='tight')
    plt.clf()

    plt.tricontourf(pts[:, 0], pts[:, 1], cb, levels=1000)
    plt.title("cmax Number")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    rect = patches.Rectangle((0.6,0.0),2.4,0.2, facecolor='white')
    ax.add_patch(rect)
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)
    plt.savefig("figures/"+name+"cmax.png", dpi=600, bbox_inches='tight')
    plt.clf()


Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/ffs05p1gridexport.h5", "output/ffs05p1solution.h5")
ffsplot(Pset, Tset, pts, vals, p, Maz, cb, "ffs05p1")

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/ffs30p1gridexport.h5", "output/ffs30p1solution.h5")
ffsplot(Pset, Tset, pts, vals, p, Maz, cb, "ffs30p1")

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/ffs05gridexport.h5", "output/ffs05solution.h5")
ffsplot(Pset, Tset, pts, vals, p, Maz, cb, "ffs05")

Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/ffs30gridexport.h5", "output/ffs30solution.h5")
ffsplot(Pset, Tset, pts, vals, p, Maz, cb, "ffs30")


