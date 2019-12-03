from espresso import espresso
from ase.build import molecule
import sys

current_pw = int(sys.argv[1])
current_k = int(sys.argv[2])

single_molecule = molecule('CO')
single_molecule.center(10)

conv_dict = {'energy':1e-6,
            'mixing':0.05,
            'mixing_mode':'local-TF',
            'maxsteps':1000,
            'diag':'cg'}

print("Calculating CO Energy")
single_molecule.calc=espresso(pw=current_pw,
                            dw=current_pw*10,
                            kpts=(1,1,1),
                            xc='PBE',
                            outdir='test_output',
                            convergence=conv_dict)

single_molecule_energy = single_molecule.get_potential_energy()

print("CO Energy: " + str(single_molecule_energy))

f = open("results.txt", "w")
f.write(str(single_molecule_energy) + '\n')
f.close()
