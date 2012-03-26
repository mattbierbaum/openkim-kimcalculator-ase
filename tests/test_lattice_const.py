import ase
from ase.structure import bulk
import numpy as np
import scipy.optimize as opt
from kimcalculator import *

symbol = 'Ar'

def energy(a, calc):
    atoms = bulk(symbol, 'fcc', a=a)
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()[0]

for m in listmodels():
    if symbol in m:
        try:
            calc = KIMCalculator(m)
            print m, " lattice const = ", opt.fmin(energy, x0=[4], args=(calc,), disp=False)
        except:
            continue
