from cosserat import CosseratRod
import numpy as np

np.set_printoptions(precision=3)

filament = CosseratRod()

ext_F=np.zeros((3,filament.N+1))
ext_F[0,filament.N] = 0.1

filament.symplectic(timespan=10,dt=0.01,method='PV',conditions=['clamp_0'],ext_F=ext_F,dissipation=1)

print(filament.x)
print('====x=====')
print(filament.Q)
print('====Q=====')
print(filament.v)
print('====v=====')
print(filament.w)
print('====w=====')
print(filament.dvdt)
print('====dvdt=====')
print(filament.dwdt)
print('====dwdt=====')
