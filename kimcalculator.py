#
#  KIM Calculator for ASE
#  
#  We make use of the SWIG python KIM interface to make use of KIM models  
#  through ASE
#  
#  Authors: Matt Bierbaum
#           YJ Chen
#
#############################################################################
#
#  Notes for work in progress:
#      1. name pointers to KIM with "km_" 
#      2. need info on neighbor list
#      3. should write a calculation_required function
#       

from kimservice import *
import kimneighborlist as kimnl
import numpy

__version__ = '0.1'
__author__  = 'Yanjiun Chen, Matt Bierbaum, Woosong Choi',

class KIMCalculator(object):

    def __init__(self,modelname):
        # correctly set up the .kim file so that is matches that
        # of modelname, i.e. will run any types
        self.testname = "test_"+modelname
        self.teststring = None

        self.modelname = modelname 
        status, km_pmdl = KIM_API_model_info(modelname)
        if KIM_STATUS_OK > status:
            KIM_API_report_error('KIM_API_model_info',status)
            raise InitializationError(self.modelname)
        KIM_API_free(km_pmdl)

        # initialize pointers for kim
        self.km_numberOfAtoms  = None
        self.km_particleCharge = None
        self.km_energy         = None
        self.km_forces         = None
        self.km_particleEnergy = None
        self.km_virial         = None
        self.km_particleVirial = None
        self.km_hessian        = None
        
        # initialize ase atoms specifications 
        self.pbc  = None 
        self.cell = None
        self.cell_orthogonal = None

        # the KIM object
        self.pkim = None

    def set_atoms(self, atoms):
        if self.pkim is not None:
            self.free_kim()
        self.init_kim(atoms)            

    def init_kim(self, atoms):
        # set the ase atoms stuff to current configuration
        self.pbc = atoms.get_pbc()
        self.cell = atoms.get_cell()
        self.cell_orthogonal = ((abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[1])) \
                               + abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[2])) \
                               + abs(numpy.dot(atoms.get_cell()[1],atoms.get_cell()[2])))<10**(-8))

        self.makeTestString(atoms)
        status, self.pkim = KIM_API_init_str(self.teststring, self.modelname)
        if KIM_STATUS_OK != status:
            KIM_API_report_error('KIM_API_init',status)
            raise InitializationError(self.modelname)
       
        natoms = atoms.get_number_of_atoms()
        ntypes = len(set(atoms.get_atomic_numbers()))

        KIM_API_allocate(self.pkim, natoms, ntypes)
        
        # set up the neighborlist as well, if necessary 
        self.uses_neighbors = uses_neighbors(self.pkim) 
        if self.uses_neighbors:
            kimnl.initialize(self.pkim)

        KIM_API_model_init(self.pkim)
       
        # get pointers to model inputs
        self.km_numberOfAtoms      = KIM_API_get_data_ulonglong(self.pkim, "numberOfParticles")
        self.km_numberAtomTypes    = KIM_API_get_data_int(self.pkim, "numberParticleTypes")
        self.km_atomTypes          = KIM_API_get_data_int(self.pkim, "particleTypes")
        self.km_coordinates        = KIM_API_get_data_double(self.pkim, "coordinates")
        if checkIndex(self.pkim,"particleCharge") >= 0:
            self.km_particleCharge = KIM_API_get_data_double(self.pkim,"particleCharge")
        if checkIndex(self.pkim, "particleSize") >= 0:
            self.km_particleSize = KIM_API_get_data_double(self.pkim,"particleSize")
    
        # check what the model calculates and get model outputs
        if checkIndex(self.pkim,"energy") >= 0:
            self.km_energy = KIM_API_get_data_double(self.pkim, "energy")
        if checkIndex(self.pkim,"forces") >= 0:
            self.km_forces = KIM_API_get_data_double(self.pkim, "forces")
        if checkIndex(self.pkim, "particleEnergy") >= 0:
            self.km_particleEnergy = KIM_API_get_data_double(self.pkim, "particleEnergy")
        if checkIndex(self.pkim, "virial") >=0:
            self.km_virial = KIM_API_get_data_double(self.pkim, "virial")
        if checkIndex(self.pkim, "particleVirial") >=0:
            self.km_particleVirial = KIM_API_get_data_double(self.pkim, "particleVirial")
        if checkIndex(self.pkim, "hessian")>=0:
            self.km_hessian = KIM_API_get_data_double(self.pkim, "hessian")
 
        
    def free_kim(self):
        if self.uses_neighbors:
            kimnl.cleanup(self.pkim)
        KIM_API_free(self.pkim)

        self.pkim = None
    
    def makeTestString(self,atoms):
        """
        makes string if it doesn't exist, if exists just keeps it as is
        """
        if self.teststring is None or self.cell_BC_changed(atoms):
            self.teststring = makekimscript(self.modelname,self.testname,atoms)

    def cell_BC_changed(self,atoms):
        """
        this function is to check whether BC has changed and cell orthogonality has changed
        because we might want to change neighbor list generator method
        """
        cell_orthogonal = ((abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[1])) \
                          + abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[2])) \
                          + abs(numpy.dot(atoms.get_cell()[1],atoms.get_cell()[2])))<10**(-8))
        if (self.pbc != atoms.get_pbc()).any() or self.cell_orthogonal != cell_orthogonal:
            return True
        else:
            return False

    def calculation_needed(self,atoms):
        """
        this checks whether or not the atoms configuration has changed and we need to recalculate..
        """
        return (self.km_energy is None or \
               (self.km_numberOfAtoms != atoms.get_number_of_atoms()) or \
               (self.km_atomTypes[:] != atoms.get_atomic_numbers()).any() or \
               (self.km_coordinates[:] != atoms.get_positions().flatten()).any() or \
               (self.pbc != atoms.get_pbc()).any() or \
               (self.cell != atoms.get_cell()).any())


    def update(self,atoms):
        """
        here we connect the KIM pointers to values in the ase atoms class
        """

        # here we only reinitialize the model if the number of Atoms / types of atoms have changed, or if the model is uninitialized
        natoms = atoms.get_number_of_atoms()
        ntypes = len(set(atoms.get_atomic_numbers()))

        if self.km_numberOfAtoms != natoms or self.km_numberAtomTypes != ntypes or self.cell_BC_changed(atoms):
            self.set_atoms(atoms)
                      
        if self.calculation_needed(atoms):
            # if the calculation is required we proceed to set the values of the standard things each model and atom class has
            self.km_numberOfAtoms[0]   = natoms 
            self.km_numberAtomTypes[0] = ntypes
            self.km_coordinates[:]     = atoms.get_positions().flatten()
            if self.km_particleCharge is not None:
                km_particleCharge[:] = atoms.get_charges()       
            
            # fill the proper chemical identifiers 
            symbols = atoms.get_chemical_symbols()
            for i in range(natoms):
                self.km_atomTypes[i] = KIM_API_get_partcl_type_code(self.pkim, symbols[i])
            
            # build the neighborlist (not a cell-based, type depends on model)
            if self.uses_neighbors:
                kimnl.set_cell(atoms.get_cell().flatten(), atoms.get_pbc().flatten().astype('int8'))
                kimnl.build_neighborlist(self.pkim)
            KIM_API_model_compute(self.pkim)
        
    def get_potential_energy(self,atoms):
        self.update(atoms)
        if self.km_energy is not None:
            return self.km_energy
        else:
            raise SupportError("energy") 
    
    def get_potential_energies(self,atoms):
        self.update(atoms)
        if self.km_particleEnergy is not None:
            particleEnergies = self.km_particleEnergy.reshape((self.km_numberOfAtoms,3))
            return particleEnergies
        else:
            raise SupportError("potential energies") 
     
    def get_forces(self,atoms):
        self.update(atoms)
        if self.km_forces is not None:
            forces = self.km_forces.reshape((self.km_numberOfAtoms,3))
            return forces
        else:
            raise SupportError("forces") 
    
    def get_stress(self,atoms):
        self.update(atoms)
        if self.km_virial is not None:
            return self.km_virial
        else:
            raise SupportError("stress")

    def get_stresses(self,atoms):
        self.update(atoms)
        if self.km_particleVirial is not None:
            return self.km_particleVirial
        else:
            raise SupportError("stress per particle")

    def get_hessian(self,atoms):
        self.update(atoms) 
        if self.km_hessian is not None:
            return self.km_hessian
        else:
            raise SupportError("hessian") 

    def __del__(self):
        """ 
        Garbage collects the KIM API objects automatically
        """
        KIM_API_free(self.pkim)


#
# KIM script maker to make a .kim-like string on the fly and 
# stores the .kim file in the appropriate directory
#
# Please refer to standard.kim for OPENKIM related standards
#

def makekimscript(modelname,testname,atoms):
    """
    input: this file takes the modelname and testname and makes an appropriate test string for the model
    and also returns it as a string
    """
    pbc = atoms.get_pbc()
    cell = atoms.get_cell()    
 
    cell_orthogonal = ((abs(numpy.dot(cell[0],cell[1])) + abs(numpy.dot(cell[0],cell[2])) + abs(numpy.dot(cell[1],cell[2])))<10**(-8))

    kimstr = "TEST_NAME :=" +testname+"\n"

    # to get model units, inputs, outputs, options we call KIM_API_model_info
    status, km_pmdl = KIM_API_model_info(modelname) 
    # this needs to be a pointer    
 
    # BASE UNIT LINES
    unit_length = KIM_API_get_unit_length(km_pmdl)
    unit_energy = KIM_API_get_unit_energy(km_pmdl)    
    unit_charge = KIM_API_get_unit_charge(km_pmdl)
    unit_temperature = KIM_API_get_unit_temperature(km_pmdl)
    unit_time = KIM_API_get_unit_time(km_pmdl)    

    kimstr += "Unit_length :=" + unit_length +"\n"
    kimstr += "Unit_energy :=" + unit_energy +"\n"
    kimstr += "Unit_charge :=" + unit_charge +"\n"
    kimstr += "Unit_temperature :=" + unit_temperature +"\n"
    kimstr += "Unit_time :=" + unit_time +"\n"


    # SUPPORTED_ATOM/PARTICLE_TYPES
    
    kimstr += "SUPPORTED_AOM/PARTICLES_TYPES: \n"    

    # check ASE atoms class for which atoms it has
    acodes = set(atoms.get_atomic_numbers())
    asymbols = set(atoms.get_chemical_symbols())

    for code, symbol in zip(list(acodes),list(asymbols)):
        kimstr += symbol +" spec "+str(code) + "\n"   

    # CONVENTIONS
    kimstr += "CONVENTIONS:\n"
    
    # note: by default the convention for python is Zero-based lists
    kimstr += "ZeroBasedLists  flag\n"
    
    # Neighbor Access methods
    # index should be non-negative and which comes first would be preferable for the model 
    index1 = checkIndex(km_pmdl, "Neigh_IterAccess")
    index2 = checkIndex(km_pmdl, "Neigh_LocaAccess")
    index3 = checkIndex(km_pmdl, "Neigh_BothAccess")

    maxindex =  max([index1, index2, index3])

    if maxindex==index1 and index1>=0:
        kimstr += "Neigh_IterAccess  flag\n"
    elif maxindex==index2 and index2>=0:
        kimstr += "Neigh_LocaAccess  flag\n"
    elif maxindex==index3 and index3>=0:
        kimstr += "Neigh_BothAccess  flag\n"
    #else:
    #   raise SupportError("Neighbor access")

    # Neighbor list and Boundary Condition (NBC) methods
    # here we can list all of the NBC methods because it will check against the model to find one that matches
    index1 = checkIndex(km_pmdl, "NEIGH_PURE_H")
    index2 = checkIndex(km_pmdl, "NEIGH_PURE_F")
    index3 = checkIndex(km_pmdl, "NEIGH_RVEC_F")
    index4 = checkIndex(km_pmdl, "MI_OPBC_H")
    index5 = checkIndex(km_pmdl, "MI_OPBC_F")
    index6 = checkIndex(km_pmdl, "CLUSTER")

    if pbc.any():
        kimstr += "NEIGH_RVEC_F flag \n"
        # we need to have RVEC if the cell is slanty
        if index3 < 0 and index4 < 0 and index5 < 0:
            raise SupportError("Periodic neighborlist")
        if cell_orthogonal:
            # we need periodicity of one variety if any pbc
            #if index3 < 0 and index4 < 0 and index5 < 0:
            #    raise SupportError("MI_OPBC_X")
            kimstr += "MI_OPBC_H flag \n"
            kimstr += "MI_OPBC_F flag \n"
        else:
            if index3 < 0:
                raise SupportError("NEIGH_RVEC_F")
    else:
        if index1 < 0 and index2 < 0 and index3 < 0 and index6 < 0:
            raise SupportError("Not periodic")
        kimstr += "NEIGH_RVEC_F flag\n"
        kimstr += "NEIGH_PURE_H flag\n"
        kimstr += "NEIGH_PURE_F flag\n"
        kimstr += "CLUSTER flag \n"        
    
    # MODEL_INPUT section
    kimstr += "MODEL_INPUT:\n"
    if checkIndex(km_pmdl,"numberOfParticles")>=0:
        kimstr +="numberOfParticles  integer  none  []\n"
    if checkIndex(km_pmdl, "numberParticleTypes")>=0:
        kimstr +="numberParticleTypes  integer  none  []\n"
    if checkIndex(km_pmdl, "particleTypes")>=0:
        kimstr +="particleTypes  integer  none  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "coordinates")>=0:
       kimstr +="coordinates  real*8  length  [numberOfParticles,3]\n" 
    if checkIndex(km_pmdl, "particleCharge")>=0:
        kimstr +="particleCharge  real*8  charge  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "particleSize")>=0:
        kimstr +="particleSize  real*8  length  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "numberContributingParticles"):
        kimstr +="numberContributingParticles  integer  none  []\n"
    
    # confused about the get_neigh and neighObject pointer.  this is test dependent.  do we include this?   
    # how do we decide whether or not to include this? Include for now. decide later
    #################################################  
    if checkIndex(km_pmdl, "get_neigh")>=0:
        kimstr += "get_neigh  method  none []\n"
    if checkIndex(km_pmdl, "neighObject")>=0:
        kimstr += "neighObject  pointer  none  []\n"
    #############################################
    if checkIndex(km_pmdl, "process_dEdr")>=0:
        kimstr +="process_dEdr  method  none  []\n"
    if checkIndex(km_pmdl, "process_d2Edr2")>=0:
        kimstr +="process_d2Edr2  method  none  []\n"
    if checkIndex(km_pmdl, "boxSideLengths")>=0:
        kimstr +="boxSideLengths  real*8  length  [3]\n"
      
    # MODEL_OUTPUT section
    kimstr += "MODEL_OUTPUT: \n"  
    if checkIndex(km_pmdl, "compute")>=0:
        kimstr += "compute  method  none  []\n"
    if checkIndex(km_pmdl, "reinit") >= 0:
        kimstr += "reinit  method  none  []\n"
    if checkIndex(km_pmdl, "destroy") >= 0:
        kimstr += "destroy  method  none  []\n"
    if checkIndex(km_pmdl, "cutoff") >= 0:
        kimstr += "cutoff  real*8  length  []\n"
    if checkIndex(km_pmdl, "energy") >= 0:
        kimstr += "energy  real*8  energy  []\n"
    if checkIndex(km_pmdl, "forces") >= 0:
        kimstr += "forces  real*8  force  [numberOfParticles,3]\n"
    if checkIndex(km_pmdl, "particleEnergy") >=0 :
        kimstr += "particleEnergy  real*8  energy  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "virial") >= 0:
        kimstr += "virial  real*8  energy  [6]\n"
    if checkIndex(km_pmdl, "particleVirial") >=0 :
        kimstr += "particleVirial  real*8  energy  [numberOfParticles,6]\n"
    if checkIndex(km_pmdl, "hessian") >= 0:
        kimstr += "hessian  real*8  pressure  [numberOfParticles,numberOfParticles,3,3]\n"
    
    return kimstr

def uses_neighbors(pkim):
    # to get model units, inputs, outputs, options we call KIM_API_model_info
    method = KIM_API_get_NBC_method(pkim) 

    if method == "CLUSTER":
       return 0
    return 1


def checkIndex(pkim,variablename):
    try:
        index = KIM_API_get_index(pkim,variablename)
    except:
        index = -1 
    return index


def listmodels():   
    import os, glob
    try:
        kimdir = os.environ['KIM_DIR']
    except:
        print "No KIM_DIR set"
        return
  
    models = [] 
    for model in glob.glob(kimdir+'MODELs/*'):
        models.append(os.path.basename(model))
    return models

class SupportError(Exception):

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)+" computation not supported by model"

class InitializationError(Exception):
            
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)+" initialization failed"
