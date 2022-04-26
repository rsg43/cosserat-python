import matplotlib.pyplot as plt
import numpy as np

class CosseratVisualisation:
	def __init__(self,filament):
		self.filament= filament
		CoM = np.mean(filament.x,axis=1)
		radius = 1.1 * (0.5 * filament.l_0 * filament.N) 
		self.box_size_x = [CoM[0]-radius,CoM[0]+radius]
		self.box_size_y = [CoM[1]-radius,CoM[1]+radius]
		self.box_size_z = [CoM[2]-radius,CoM[2]+radius]

	def simple_lineplot(self):
		x = self.filament.x[0,:]
		y = self.filament.x[1,:]
		z = self.filament.x[2,:]
		fig = plt.figure
		ax = plt.axes(projection='3d')
		ax.plot3D(x,y,z)
		ax.set_xlim(self.box_size_x)
		ax.set_ylim(self.box_size_y)
		ax.set_zlim(self.box_size_z)
		plt.show()

	def cylinder_coords(self,radius,length,orientation,CoM):
		z = np.linspace(-0.5 * length, 0.5 * length, 100)
		theta = np.linspace(0, 2*np.pi, 100)
		theta_grid, z_grid=np.meshgrid(theta, z)
		x_grid = radius*np.cos(theta_grid)
		y_grid = radius*np.sin(theta_grid)
		x_grid += CoM[0]
		y_grid += CoM[1]
		z_grid += CoM[2]
		return x_grid, y_grid, z_grid

	def cylinder_plot(self):
		radius = np.sqrt(self.filament.A)
		length = 1

		x_grid, y_grid, z_grid = cylinder_coords(radius,length,orientation,CoM)


