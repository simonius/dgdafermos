# This module contains all needed functionality to import the exported solutions from HDF5
# file to produce images using the matplotlib. 

import numpy as np
import h5py
import matplotlib.pyplot as plt
import math
from scipy.interpolate import griddata


from euler import *

Psetname = "Points"
Tsetname = "Triangles"
Asetname = "Adjacents"
CPsetname = "ColPoints"

# Imports a solution
# gfname is the name of the grid export HDF 5 file
# solfname is the name of the solution export HDF 5 file
def SolImport(gfname, solfname):

    # Loading of the grid from the grid export file
    f = h5py.File(gfname, "r")

    # The points in the grid 
    Pset = np.transpose(f[Psetname])
    # The triangles - an array pointing into Pset
    Tset = np.transpose(f[Tsetname]) - 1 # Convert to Zero based indexes!
    # The adjacencies
    Aset = np.transpose(f[Asetname]) - 1 # Convert to zero based indexes!
    # The collocation points in the reference triangle
    CPset = np.transpose(f[CPsetname])

    # loading the solution from the solution export
    f = h5py.File(solfname, "r")
    u = np.transpose(f["u"])

    # buildup of the grid points and solution 
    Ntr = Tset.shape[0]
    Ncp = CPset.shape[0]
    pts = np.zeros([Ntr*Ncp, 2])
    PTpts = np.zeros([Ntr, Ncp, 2])

    rads = np.zeros([Ntr*Ncp])
    vals = np.zeros([Ntr*Ncp, 4])
    p = np.zeros([Ntr*Ncp])
    print("Using ", Ntr, "triangles")
    print("File of nodal values has ", u.shape[0], " triangles")
    print("We are using ", Ncp, "collocation points per triangle")
    for tr in range(Ntr):
        x1 = Pset[Tset[tr, 0], :]
        x2 = Pset[Tset[tr, 1], :]
        x3 = Pset[Tset[tr, 2], :]
        ToLK = np.zeros([2, 2])
        ToLK[:, 0] = x2-x1
        ToLK[:, 1] = x3-x1
        for cp in range(Ncp):
            pts[cp + Ncp*tr, :] = Pset[Tset[tr, 0], :] + np.matmul(ToLK,CPset[cp, :])
            PTpts[tr, cp, :] = Pset[Tset[tr, 0], :] + np.matmul(ToLK,CPset[cp, :])
            rads[cp + Ncp*tr] = np.linalg.norm(pts[cp + Ncp*tr, :])
            for k in range(4):
                if np.isnan(u[tr, cp, k]):
                    vals[cp + Ncp*tr, k] = 0.0
                else:
                    vals[cp + Ncp*tr, k] = u[tr, cp, k]
            p[cp + Ncp*tr] = press(vals[cp + Ncp*tr, :])


    return Pset, Tset, pts, vals, p, PTpts, u

# Plots a grid saved in Pset and Tset
def PlotGrid(Pset, Tset):
    for tr in range(Tset.shape[0]):
        kords = np.zeros([4, 2])
        for k in range(3):
            kords[k, :] = Pset[Tset[tr, k], :]
        kords[3, :] = Pset[Tset[tr, 0], :]
        plt.plot(kords[:, 0], kords[:, 1], color="k", linewidth=0.3)


