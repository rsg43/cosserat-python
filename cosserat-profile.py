from cosserat import CosseratRod
import cProfile

filament = CosseratRod()
cProfile.run('filament.symplectic()')