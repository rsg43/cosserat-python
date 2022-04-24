from cosserat import CosseratRod
from cosserat_visualisation import CosseratVisualisation
import numpy as np

#set precision for print outputs fr numpy arrays
np.set_printoptions(precision=3)
#intialise Cosserat filament 
filament = CosseratRod()
plotter = CosseratVisualisation(filament)
#preallocate external force array
ext_F=np.zeros((3,filament.N+1))
#set force at end of filament 
ext_F[0,filament.N] = 2
#simulate filament using sympletic, with clamped end and dissipation
filament.symplectic(timespan=100,conditions=['clamp_0'],ext_F=ext_F,dissipation=1)
#print position of filament after simulation
print(filament.x)
print('====x=====')

plotter.simple_lineplot()