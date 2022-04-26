from cosserat import CosseratRod
import numpy as np

class CosseratMultipleFilament:
	def __init__(self, filaments):
		self.filaments = filaments
		self.num_filaments = len(filaments)

	def linkers(self):
		pass

	def symplectic(self,timespan=10,dt=0.01,method='PEFRL',dissipation=0,conditions=[]):
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

		#Set number of steps in simulation
		numsteps = int(timespan / dt)
		#Run simulation loop
		for ii in range(numsteps):
            #Loop through integrator coefficents
			for jj in range(max(len(a),len(b))):
                #check to see if there is a first velocity step
				if a[jj] != 0:
					self.linkers()
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


filaments = [None] * 10
for ii in range(10):
	filaments[ii] = CosseratRod()


filament_store = CosseratMultipleFilament(filaments)
filament_store.symplectic()