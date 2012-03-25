import ase
import numpy as np
import scipy.optimize as opt
from kimcalculator import *

calc = KIMCalculator("ex_model_Ar_P_LJ")
a = 5
def energy(a):
    atoms = ase.Atoms( symbols='Ar', positions=np.array([[0,0,0]]), pbc=[(1,1,1)], cell=np.array([[a+1e-6,0,0],[0,a,0],[0,0,a]]) )
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()[0]
    
min = opt.fmin(energy, [1])
print min
