from cosserat import CosseratRod
import cProfile
import numpy as np

filament = CosseratRod()
ext_F=np.zeros((3,filament.N+1))
ext_F[0,filament.N] = 2

cProfile.run('filament.symplectic(timespan=10,conditions=[\'clamp_0\'],ext_F=ext_F,dissipation=1)')

def logm_test():
	for i in range(100000):
		filament.logm(filament.Q[:,:,3])

cProfile.run('logm_test()')

def expm_test():
	for i in range(100000):
		filament.expm(filament.w[:,3], 0.01)

cProfile.run('expm_test()')
