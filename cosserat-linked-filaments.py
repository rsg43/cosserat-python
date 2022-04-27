from cosserat import CosseratRod
from cosserat_multiple_filaments import CosseratMultipleFilament
from cosserat_visualisation import CosseratVisualisation


num_fils = 3
filaments = [None] * num_fils
for ii in range(num_fils):
	filaments[ii] = CosseratRod()
links = [[[0,0],[1,3]],[[1,7],[2,5]],[[2,3],[0,4]],[[1,9],[0,8]],[[2,7],[0,5]]]


filament_store = CosseratMultipleFilament(filaments)
filament_store.symplectic(timespan=10,links=links,dissipation=1)


print(filaments[1].x)

plotter = CosseratVisualisation(filaments[0])
plotter.simple_lineplot()
