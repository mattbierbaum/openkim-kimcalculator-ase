from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
from kimcalculator import *
from numpy import *

ar = FaceCenteredCubic('Ar', pbc=[(1,1,1)], directions=[[1,0,0],[0,1,0],[1,1,1]], size=[4,4,4])

for m in listmodels():
    if "Ar" in m:
        try:
            print "Try model: ", m
            ar.set_pbc([(1,1,1)])
            calc1 = KIMCalculator(m)
            ar.set_calculator(calc1)
           
            N = ar.get_number_of_atoms()
            bulk_energy = ar.get_potential_energy()[0] / N
            
            ar.set_pbc([(1,1,0)])
            surf_energy = ar.get_potential_energy()[0] / N
            print "\tsurface = ", surf_energy - bulk_energy  
            print      
        except SupportError, e:
            print "\tskipping ", m
            print "\t", e, " ...\n"
            continue 
