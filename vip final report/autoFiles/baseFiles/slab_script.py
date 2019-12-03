from espresso import espresso
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import surface
from ase.spacegroup import crystal
from ase.io import read
import sys

current_pw = int(sys.argv[1])
current_k = int(sys.argv[2])

a=3.923 #angstrom
Pt=crystal(['Pt'],basis=[(0,0,0)],spacegroup=225,cellpar=[a,a,a,90,90,90])

single_slab=surface(Pt,(1,1,1),2)
single_slab.center(vacuum=10,axis=2)
single_slab.set_calculator(EMT())
constraint = FixAtoms(mask=[False, True, True, False, False, True, True, False])
single_slab.set_constraint(constraint)

dyn = QuasiNewton(single_slab, trajectory='modified_slab.traj')
dyn.run(fmax=1)

modified_single_slab = read('modified_slab.traj')

conv_dict = {'energy':1e-6,
            'mixing':0.05,
            'mixing_mode':'local-TF',
            'maxsteps':1000,
            'diag':'cg'}

print("Calculating slab Energy")

modified_single_slab.calc=espresso(pw=current_pw,
                       dw=current_pw*10,
                        kpts=(current_k,current_k,1),
                       xc='PBE',
                       outdir='test_output',
                       convergence=conv_dict)

single_slab_energy = modified_single_slab.get_potential_energy()

print("Slab Energy: " + str(single_slab_energy))

f = open("results.txt", "w")
f.write(str(single_slab_energy) + '\n')
f.close()
