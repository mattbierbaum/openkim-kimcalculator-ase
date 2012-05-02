from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
from kimcalculator import *
from numpy import *

N = 2
ar = FaceCenteredCubic('Ar', pbc=[(1,1,1)], directions=[[1,0,0],[0,1,0],[0,0,1]], size=[N,N,N])
print ar.get_cell()
#view(ar) 

calc1 = KIMCalculator("ex_model_Ar_P_LJ")
ar.set_calculator(calc1)
kim_energy = ar.get_potential_energy()
print "kim energy = ", kim_energy 
