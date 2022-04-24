import matplotlib.pyplot as plt

class CosseratVisualisation:
	def __init__(self,filament):
		self.filament=filament

	def simple_lineplot(self):
		x = self.filament.x[0,:]
		y = self.filament.x[1,:]
		z = self.filament.x[2,:]
		fig = plt.figure
		ax = plt.axes(projection='3d')
		ax.plot3D(x,y,z)
		plt.show()

	def cylinder_coords(self,radius,length,orientation):
		# z = np.linspace(0, height_z, 50)
	 #    theta = np.linspace(0, 2*np.pi, 50)
	 #    theta_grid, z_grid=np.meshgrid(theta, z)
	 #    x_grid = radius*np.cos(theta_grid) + center_x
	 #    y_grid = radius*np.sin(theta_grid) + center_y
	 #    return x_grid,y_grid,z_grid
		pass

	def cylinder_plot(self):
		pass




