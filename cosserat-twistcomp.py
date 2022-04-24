from cosserat import CosseratRod
import numpy as np

#set precision for print outputs fr numpy arrays
np.set_printoptions(precision=3)
#intialise Cosserat filament 
filament = CosseratRod()
conditions = ['comp_clamp_0','comp_clamp_N','twist_clamp_0','twist_clamp_N','noise_no_ext']
#simulate filament using sympletic, with clamped end and dissipation
filament.symplectic(timespan=100,conditions=conditions,dissipation=1)
#print position of filament after simulation
print(filament.x)
print('====x=====')
