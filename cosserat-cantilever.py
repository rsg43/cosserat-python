from cosserat import CosseratRod
import numpy as np

np.set_printoptions(precision=3)

filament = CosseratRod()

ext_F=np.zeros((3,filament.N+1))
ext_F[0,filament.N] = 0.1

filament.symplectic(timespan=0.5,method='VV',conditions=['clamp_0'],ext_F=ext_F)
print(filament.x)
print(filament.v)
print(filament.w)