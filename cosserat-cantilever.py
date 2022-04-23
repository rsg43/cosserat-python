from cosserat import CosseratRod
import numpy as np

np.set_printoptions(precision=3)

filament = CosseratRod()

ext_F=np.zeros((3,filament.N+1))
ext_F[0,filament.N] = 0.01

print(filament.x)
print(filament.sigma)
print(filament.kappa)
filament.symplectic(timespan=10,method='PEFRL',conditions=['clamp_0'],ext_F=ext_F,dissipation=0.001)
print(filament.x)
print(filament.sigma)
print(filament.kappa)
print(filament.Q)

