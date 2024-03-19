# This python module allows to calculate physical variables,
# sound speeds, mach numbers and speed bounds from conserved variables. 

import numpy as np
import math

# Calculates the pressure
def press(u):
    gamma = 1.4
    try:
        p = (gamma-1)*(u[3]-0.5*(u[1]**2 + u[2]**2)/u[0])
    except:
        p = 0.0
    return p

# calculates the sound speed
def a(u):
    gamma = 1.4
    try:
        a = math.sqrt(gamma*press(u)/u[0])
    except:
        a = 0.0
    return a

# calculates the Mach number
def Ma(u):
    try:
        v = u[1:3]/u[0]
        Ma = np.linalg.norm(v)/a(u)
    except:
        Ma = 0.0
    return Ma

# Calculates a speed bound
def cbound(u):
    try:
        v = u[1:3]/u[0]
        va = np.linalg.norm(v) + a(u)
    except:
        va = 0.0
    return va
