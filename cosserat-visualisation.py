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


