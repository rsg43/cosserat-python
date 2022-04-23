from cosserat import CosseratRod

filament = CosseratRod()
print(filament.x)
print(filament.Q)
print(filament.B)
print(filament.x[:,3])
print(filament.expm(filament.x[:,0]))
filament.update_Q(0.001)
print(filament.logm(filament.Q[:,:,3]))
print(filament.sigma)
filament.update_sigma()
print(filament.sigma)
filament.update_kappa()
print(filament.update_acceleration())
filament.symplectic()
print(filament.x)