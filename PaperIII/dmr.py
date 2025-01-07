# This file produces the plots of the forward facing step problem
# from the saved HDF5 files in output

from solimport import *
from matplotlib import patches
from schlieren import *

import matplotlib.tri as tri



def ffsplot(Pset, Tset, pts, vals, p, ptpts, u,  name):
    # Calculation of the schlieren map
    smap = schlieren.ptschl(ptpts, u, 1)

    gridp = True
    pcontourf = False
    rcontourf = True
   
    if gridp:
        PlotGrid(Pset, Tset)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel("x")
        plt.ylabel("y")

        plt.savefig("figures/"+name+"grid.png", dpi=600, bbox_inches='tight')
        plt.clf()

    if rcontourf:
        rs = np.sort(vals[:, 0])
        ci = int(math.floor(math.log(len(rs))))
        minval = rs[ci]
        maxval = rs[-ci]
        cl = np.linspace(minval, maxval, num=1000)
        plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=cl)
        plt.title("Density")
        plt.xlabel("x")
        plt.ylabel("y")
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)

        plt.savefig("figures/"+name+"rho.png", dpi=600, bbox_inches='tight')
        plt.clf()


#   Schlieren plot
    plt.tricontourf(pts[:, 0], pts[:, 1], np.sqrt(smap), levels=1000, cmap='Greys')
    plt.title("Density, Schlieren image")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig("figures/"+name+"smap.png", dpi=600, bbox_inches='tight')
    plt.clf()
    
    if pcontourf:
        plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
        plt.title("Pressure")
        plt.xlabel("x")
        plt.ylabel("y")
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)
        plt.savefig("figures/"+name+"p.png", dpi=600, bbox_inches='tight')
        plt.clf()


Pset, Tset, pts, vals, p, ptpts, u  = SolImport("output/dmrsolp1gridexport.h5", "output/dmrsolp1solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "dmrsolp1")

Pset, Tset, pts, vals, p, ptpts, u  = SolImport("output/dmrsolp2gridexport.h5", "output/dmrsolp2solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "dmrsolp2")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dmrsolp3gridexport.h5", "output/dmrsolp3solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u,  "dmrsolp3")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dmrsolp5gridexport.h5", "output/dmrsolp5solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "dmrsolp5")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dmrsolp7gridexport.h5", "output/dmrsolp7solution.h5")
ffsplot(Pset, Tset, pts, vals, p, ptpts, u, "dmrsolp7")
