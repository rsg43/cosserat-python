from cosserat import CosseratRod
import numpy as np

np.set_printoptions(precision=3)

filament = CosseratRod()

ext_F=np.zeros((3,filament.N+1))
ext_F[0,filament.N] = 2
ext_F[0,filament.N//2] = -2
ext_F[0,0] = 2

filament.symplectic(timespan=10,ext_F=ext_F,dissipation=1)

print(filament.x)
print('====x=====')