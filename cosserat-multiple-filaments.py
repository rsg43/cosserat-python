from cosserat import CosseratRod
from cosserat_visualisation import CosseratVisualisation
import numpy as np

#set precision for print outputs fr numpy arrays
np.set_printoptions(precision=3)
#intialise Cosserat filament 
filament = CosseratRod()
#initialise plotter class
plotter = CosseratVisualisation(filament)
#simulate filament using sympletic, with clamped end and dissipation
filament.symplectic(timespan=100,dissipation=1)
#print position of filament after simulation
print(filament.x)
print('====x=====')
#make simple lineplot of filament centreline
plotter.simple_lineplot()