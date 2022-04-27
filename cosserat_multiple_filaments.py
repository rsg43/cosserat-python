from cosserat import CosseratRod
import numpy as np

class CosseratMultipleFilament:
	def __init__(self, filaments):
		self.filaments = filaments
		self.num_filaments = len(filaments)

	def linkers_setup(self,links):
		self.linker_length = 5
		self.linker_stiffness = 100
		self.links = links

		for link in links:
			if link[0][0] >= self.num_filaments \
			and link[1][0] >= self.num_filaments \
			and link[0][1] >= self.filaments[link[0][0]].N \
			and link[1][1] >= self.filaments[link[1][0]].N:
				raise ValueError(f'{i} is an invalid bond')

	def linkers_force(self):
		for nn in range(self.num_filaments):
			force_scale = 0.01
			self.filaments[nn].ext_F = np.random.normal(size=(3,self.filaments[nn].N+1)) * force_scale
		for link in self.links:
			x0 = np.copy(self.filaments[link[0][0]].x[:,link[0][1]])
			x1 = np.copy(self.filaments[link[1][0]].x[:,link[1][1]])
			f = self.linker_stiffness * (x1 - x0) * (1 - self.linker_length / np.sqrt(np.einsum('i,i->',x1 - x0,x1 - x0)))
			self.filaments[link[0][0]].ext_F[:,link[0][1]] = f
			self.filaments[link[1][0]].ext_F[:,link[1][1]] = -f

	def symplectic(self,timespan=50,dt=0.01,method='PEFRL',dissipation=0,conditions=[],links=None):
		#Set integrating method as position extended Forrest-Ruth like   
		if method == 'PEFRL':
		    xi = 0.1786178958448091
		    lmbda = -0.2123418310626054
		    chi = -0.06626458266981849
		    a = np.array([xi,chi,1-2*(xi+chi),chi,xi])
		    b = np.array([0.5*(1-2*lmbda),lmbda,lmbda,0.5*(1-2*lmbda)])
		#Set integrating method as velocity Verlet
		elif method == 'VV':
		    a = np.array([0.5,0.5])
		    b = np.array([1])
		#Set integrating method as position Verlet
		elif method == 'PV':
		    a = np.array([0,1])
		    b = np.array([0.5,0.5])
		#Raise error if method not recognised
		else:
		    raise ValueError('Incompatible integrator specified')
		
		if isinstance(links, list):
			self.linkers_setup(links)
		#Set number of steps in simulation
		numsteps = int(timespan / dt)
		#Run simulation loop
		for ii in range(numsteps):
            #Loop through integrator coefficents
			for jj in range(max(len(a),len(b))):
                #check to see if there is a first velocity step
				if a[jj] != 0:
					if isinstance(links, list):
						self.linkers_force()
					for nn in range(self.num_filaments):
	                    #update accelerations
						self.filaments[nn].update_acceleration()
	                    #update velocities
						self.filaments[nn].update_v(a[jj] * dt, dissipation)
						self.filaments[nn].update_w(a[jj] * dt, dissipation)
	                    #enforce boundary conditions
						self.filaments[nn].update_conditions(conditions,numsteps,ii,dt)
                #check to see if there is a final position step
				if jj < len(b):
					for nn in range(self.num_filaments):
	                    #update positions and orientations
						self.filaments[nn].update_x(b[jj] * dt)
						self.filaments[nn].update_Q(b[jj] * dt)

