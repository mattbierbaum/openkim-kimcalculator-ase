from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
from kimcalculator import *
from numpy import *

symbol = 'Ar'
ar = FaceCenteredCubic(symbol, pbc=[(1,0,0)], directions=[[1,0,0],[0,1,0],[0,0,1]], size=[5,5,5])

for m in listmodels():
    if symbol in m:
        try:
            print "Try model: ", m
            ar.set_pbc([(1,1,1)])
            calc1 = KIMCalculator(m)
            ar.set_calculator(calc1)
           
            N = ar.get_number_of_atoms()
            ar.set_pbc([(1,1,1)])
            bulk_energy = ar.get_potential_energy()[0] / N

            ar.set_pbc([(1,1,0)])
            surf_energy = ar.get_potential_energy()[0] / N
            print "\tsurface = ", surf_energy - bulk_energy  

            #virial = ar.get_stresses()
            #print "\tvirial = ", virial
            print      
        except SupportError, e:
            print "\tskipping ", m, "\n\t", e, " ...\n"
            continue 


