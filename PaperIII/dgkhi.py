# This file produces the plots of the Kelvin-Helmholtz problem
# from the saved HDF5 files in output

from solimport import *
from matplotlib import patches


def khiplot(Pset, Tset, pts, vals, p, ptpts, u, name):
  
    # Plots the grid
    PlotGrid(Pset, Tset)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("figures/"+name+"grid.png", dpi=600, bbox_inches='tight')
    plt.clf()

    # Plot the density
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



    # plot the magnified versions of the solution with a grid overlay
    PlotGrid(Pset, Tset)
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)

    plt.xlim(0.5, 1.0)
    plt.ylim(0.5, 1.0)
    plt.savefig("figures/"+name+"rhocloseL.png", dpi=600, bbox_inches='tight')
    plt.xlim(3.0, 3.5)
    plt.ylim(0.5, 1.0)
    plt.savefig("figures/"+name+"rhocloseR.png", dpi=600, bbox_inches='tight')
    plt.clf()
    
    # Plot the velocity in x
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 1] / vals[:, 0], levels=1000)
    plt.title("Velocity vx")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)

    plt.savefig("figures/"+name+"vx.png", dpi=600, bbox_inches='tight')
    plt.clf()

    # Plot the velocity in y
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 2] / vals[:, 0], levels=1000)
    plt.title("Velocity vy")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)

    plt.savefig("figures/"+name+"vy.png", dpi=600, bbox_inches='tight')
    plt.clf()


    # Plot the pressure
    plt.tricontourf(pts[:, 0], pts[:, 1], p, levels=1000)
    plt.title("Pressure")
    plt.xlabel("x")
    plt.ylabel("y")
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(fraction=0.046, pad=0.04, shrink=0.5)
    plt.savefig("figures/"+name+"p.png", dpi=600, bbox_inches='tight')
    plt.clf()


Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dgkhip1gridexport.h5", "output/dgkhip1solution.h5")
khiplot(Pset, Tset, pts, vals, p, ptpts, u, "dgkhip1")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dgkhip2gridexport.h5", "output/dgkhip2solution.h5")
khiplot(Pset, Tset, pts, vals, p, ptpts, u, "dgkhip2")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dgkhip3gridexport.h5", "output/dgkhip3solution.h5")
khiplot(Pset, Tset, pts, vals, p, ptpts, u, "dgkhip3")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dgkhip5gridexport.h5", "output/dgkhip5solution.h5")
khiplot(Pset, Tset, pts, vals, p, ptpts, u, "dgkhip5")

Pset, Tset, pts, vals, p, ptpts, u = SolImport("output/dgkhip7gridexport.h5", "output/dgkhip7solution.h5")
khiplot(Pset, Tset, pts, vals, p, ptpts, u, "dgkhip7")





