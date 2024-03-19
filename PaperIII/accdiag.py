# This file is used to look at the errors produced by the accuracy
# test

from solimport import *

def ErrorPlot(name, Pset, Tset, pts, vals, p, Maz):
    plt.tricontourf(pts[:, 0], pts[:, 1], vals[:, 0], levels=1000)
    plt.colorbar()
    plt.title("Density Error")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.savefig("figures/"+name+"rho.png", dpi=600,bbox_inches='tight')
    plt.clf()


Pset, Tset, pts, vals, p, Maz, cb = SolImport("output/acc1gridexport.h5", "output/acc1solution.h5")
ErrorPlot("Acc1", Pset, Tset, pts, vals, p, Maz)
