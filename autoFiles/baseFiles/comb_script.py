from espresso import espresso
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import molecule, surface, add_adsorbate
from ase.spacegroup import crystal
from ase.io import read, write
import sys

current_pw = int(sys.argv[1])
current_k = int(sys.argv[2])

a=3.923 #angstrom
Pt=crystal(['Pt'],basis=[(0,0,0)],spacegroup=225,cellpar=[a,a,a,90,90,90])

h = 1.85

conv_dict = {'energy':1e-6,
            'mixing':0.05,
            'mixing_mode':'local-TF',
            'maxsteps':1000,
            'diag':'cg'}

Pt_slab=surface(Pt,(1,1,1),2)
Pt_slab.center(vacuum=10,axis=2)
CO_molecule = molecule('CO')
CO_molecule.center(10)

CO_molecule.set_calculator(EMT())
Pt_slab.set_calculator(EMT())

add_adsorbate(Pt_slab, CO_molecule, h, position=(1, 1))

constraint = FixAtoms(mask=[False, True, True, False, False, True, True, False])
Pt_slab.set_constraint(constraint)

dyn = QuasiNewton(Pt_slab, trajectory='Pt+CO.traj')
dyn.run(fmax=3)

dft_slab = read('Pt+CO.traj')

dft_slab.calc=espresso(pw=current_pw,
                       dw=current_pw*10,
                       kpts=(current_k,current_k,1),
                       xc='PBE',
                       outdir='test_output',
                       convergence=conv_dict)

combined_energy = dft_slab.get_potential_energy()

f = open("results.txt", "w")
f.write(str(combined_energy) + '\n')
f.close()
