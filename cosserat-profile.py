from cosserat import CosseratRod
import cProfile

filament = CosseratRod()
cProfile.run('filament.symplectic(timespan=100)')

def logm_test():
	for i in range(1000000):
		filament.logm(filament.Q[:,:,3])


cProfile.run('logm_test()')
