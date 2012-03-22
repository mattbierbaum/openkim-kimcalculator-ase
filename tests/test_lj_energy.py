from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
from kimcalculator import *
from numpy import *

cutoff = 8.15
epsilon = 0.0104
sigma = 3.4

import numpy as np
class LennardJones:
    def __init__(self, epsilon=0.0104, sigma=3.4, cutoff=8.15):
        self.epsilon = epsilon
        self.sigma = sigma
        self.cutoff = cutoff
        self.positions = None
        
        c6 = (self.sigma**2 / self.cutoff**2)**3
        c12 = c6**2
        self.shift = 4 * self.epsilon*(c12-c6)

    def update(self, atoms):
        assert not atoms.get_pbc().any()
        if (self.positions is None or
            (self.positions != atoms.get_positions()).any()):
            self.calculate(atoms)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self._forces

    def get_stress(self, atoms):
        return np.zeros((3, 3))
    
    def calculate(self, atoms):
        positions = atoms.get_positions()
        self.energy = 0.0
        self._forces = np.zeros((len(atoms), 3))
        for i1, p1 in enumerate(positions):
            for i2, p2 in enumerate(positions[:i1]):
                diff = p2 - p1
                d2 = np.dot(diff, diff)
                c6 = (self.sigma**2 / d2)**3
                c12 = c6**2
                if d2 < self.cutoff**2:
                    self.energy += 4 * self.epsilon * (c12 - c6) - self.shift
                F = 24 * self.epsilon * (2 * c12 - c6) / d2 * diff
                self._forces[i1] -= F
                self._forces[i2] += F
        self.positions = positions.copy()



N = 4
ar = FaceCenteredCubic('Ar', pbc=[(0,0,0)], directions=[[1,8,3],[0,1,0],[1,1,1]], size=[N,N,N])
print ar.get_cell()
#view(ar) 

calc1 = KIMCalculator("ex_model_Ar_P_LJ")
ar.set_calculator(calc1)
kim_energy = ar.get_potential_energy()[0]
print "kim energy = ", kim_energy 

calc2 = LennardJones(epsilon=epsilon, sigma=sigma, cutoff=cutoff)
ar.set_calculator(calc2)
ase_energy = ar.get_potential_energy()
print "ase energy = ", ase_energy 

print "difference = ", kim_energy - ase_energy
