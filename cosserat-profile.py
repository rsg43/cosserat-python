from cosserat import CosseratRod
import cProfile

filament = CosseratRod()
cProfile.run('filament.symplectic(timespan=10,conditions=[\'clamp_0\'])')

# def logm_test():
# 	for i in range(100000):
# 		filament.logm(filament.Q[:,:,3])


# cProfile.run('logm_test()')
