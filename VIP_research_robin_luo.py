#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:54:31 2019

@author: robin
"""

from ase.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.build import molecule
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import surface, fcc111, add_adsorbate
from ase.io import read
from ase.optimize import QuasiNewton
from gpaw import GPAW, PW




# instruction 1 below
#build the Pt crystal structure
a=3.923
Pt=crystal(['Pt'],basis=[(0,0,0)],spacegroup=225,
           cellpar=[a,a,a,90,90,90])

#build the (111) surface slab
Pt_111=surface(Pt,(1,1,1),2)
Pt_111.center(vacuum=10,axis=2)

#repeat slab
Pt_111_repeat=Pt_111.repeat((1,2,1))


#save .traj file
write('Pt_111.traj',Pt_111_repeat)
write('POSCAR',Pt_111) #in case the .traj file cannot be read
view(Pt_111_repeat)



#instruction 2 below

atoms=molecule('CO2')
atoms.center(10)

#write the .traj file
write('ads_CO.traj',atoms)
write('POSCAR',atoms)


#instruction 3 below

#the paths to your .traj file
Pt_slab_path = 'Pt_111.traj'
ads_CO_path = 'ads_CO.traj'

#the height in angstroms the gas will be above the slab
h = 4.83

#read in the .traj files and put them into variables
Pt_slab = read(Pt_slab_path)
CO_molecule = read(ads_CO_path)

#these following lines will perform a calculation to slightly optimize the #position of the gas on the slab, but is no substitute for DFT calculations
CO_molecule.set_calculator(EMT())
Pt_slab.set_calculator(EMT())
E_Pt=Pt_slab.get_total_energy()

add_adsorbate(Pt_slab,CO_molecule,h,position=(2,1),mol_index=1)
#run the EMT optimization and write to a file


dyn = QuasiNewton(Pt_slab, trajectory='Pt+CO.traj')
dyn.run(fmax=0.15)
E_CO=CO_molecule.get_total_energy()
E_total=Pt_slab.get_total_energy()


#open a GUI to inspect the adsorbed gas
view(Pt_slab)
print("the Platium energy before:")
print(E_Pt)
print("the CO energy before:")
print(E_CO)
print("the energy after:")
print(E_total)

#instruction 4 below
slab_ads=read('Pt+CO.traj')
slab_ads.calc=GPAW(mode=PW(300),
                       kpts=(1,1,1),
                       xc='PBE',
                  spinpol=True,
                  convergence={'energy':1e-6})
relax_slab_ads=QuasiNewton(slab_ads,
                           logfile='opt.log',
                           trajectory='opt.traj',
                           restart='opt.pckl')
relax_slab_ads.run(fmax=0.15)
