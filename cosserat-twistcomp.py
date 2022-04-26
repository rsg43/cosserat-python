from cosserat import CosseratRod
from cosserat_visualisation import CosseratVisualisation
import numpy as np
import matplotlib.pyplot as plt

#set precision for print outputs fr numpy arrays
np.set_printoptions(precision=3)
#intialise Cosserat filament 
filament = CosseratRod()
#initialise plotter class
plotter = CosseratVisualisation(filament)
#set conditions for simulations
conditions = ['comp_clamp_0','comp_clamp_N','twist_clamp_0','twist_clamp_N','noise_no_ext']
#simulate filament using sympletic, with clamped end and dissipation
filament.symplectic(timespan=1,conditions=conditions,dissipation=1)
#print position of filament after simulation
print(filament.x)
print('=====x=====')
#make simple lineplot of filament centreline
plotter.simple_lineplot()

#set params for cylinder test
radius,length,orientation,CoM = 5,10,None,[0,0,0]
#reset axes
ax = plt.axes(projection='3d')
#create surface plot
x_grid,y_grid,z_grid = plotter.cylinder_coords(radius,length,orientation,CoM)
ax.plot_surface(x_grid,y_grid,z_grid)
#show figure
plt.show()
