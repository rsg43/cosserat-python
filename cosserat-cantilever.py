from cosserat import CosseratRod

filament = CosseratRod()
filament.symplectic(timespan=10,conditions=['clamp_0'])
print(filament.x)