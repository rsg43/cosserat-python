from cosserat import CosseratRod
import numpy as np

#set precision for print outputs fr numpy arrays
np.set_printoptions(precision=3)
#intialise Cosserat filament
filament = CosseratRod()
#preallocate external force array
ext_F=np.zeros((3,filament.N+1))
#set (balanced) forces acting on ends and middle of filament
ext_F[0,filament.N] = 1
ext_F[0,filament.N//2] = -2
ext_F[0,0] = 1
#simulate filament using sympletic, with dissipation
filament.symplectic(timespan=100,ext_F=ext_F,dissipation=1)
#print position of filament after simulation
print(filament.x)
print('====x=====')