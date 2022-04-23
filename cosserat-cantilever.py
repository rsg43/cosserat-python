from cosserat import CosseratRod
import numpy as np

np.set_printoptions(precision=3)

filament = CosseratRod()

ext_F=np.zeros((3,filament.N+1))
ext_F[0,filament.N] = 1

filament.symplectic(timespan=100,dt=0.001,method='PV',conditions=['clamp_0'],ext_F=ext_F,dissipation=10)


