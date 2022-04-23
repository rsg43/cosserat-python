import numpy as np

class CosseratRod:

    def __init__(self, segments=10, seg_length=10.0, position=[[0.0],[0.0],[0.0]]):
        #Basic geometric properties
        self.N = segments
        self.l_0 = seg_length
        self.A = 38.5
        self.rho = 0.01

        #Current coordinates, velocities, strains, Voronoi, extensions
        self.x = np.zeros((3,self.N+1))
        self.x[2,:] = np.array([i*self.l_0 for i in range(self.N+1)])
        self.x += np.array(position)
        self.v = np.zeros((3,self.N+1))
        self.Q = np.zeros((3,3,self.N))
        for ii in range(self.N):
            self.Q[:,:,ii] = np.identity(3)
        self.w = np.zeros((3,self.N))
        self.kappa = np.zeros((3,self.N-1))
        self.sigma = np.zeros((3,self.N))
        self.l = self.x[:,1:] - self.x[:,:self.N]
        self.l_norm = np.linalg.norm(self.l,axis=0)
        self.D_0 = 0.5 * (self.l_norm[1:] + self.l_norm[:self.N-1])
        self.D = np.copy(self.D_0)
        self.e = self.l_norm / self.l_0
        self.ee = self.D / self.D_0

        #Mass, inertia, rigidities
        self.m = (seg_length * self.A * self.rho) * np.ones((1,self.N+1))
        self.m[0,0] *= 0.5
        self.m[0,self.N] *= 0.5
        self.I = ((self.m[0,1] * self.A ** 2) / (4 * np.pi)) * np.diag(np.array([1.0,1.0,2.0]))
        self.B = np.diag(np.array([39000.0, 39000.0, 13000.0]))
        self.S = np.diag(np.array([16000.0, 16000.0, 34000.0]))
        
    def expm(self,u):
        theta = np.linalg.norm(u)
        if theta > 0.0:
            u /= theta
            U = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
            return np.identity(3) + np.sin(theta) * U + np.cos(theta) * U * U
        return np.identity(3)
    
    def logm(self,R):
        theta = np.arccos(0.5 * (np.einsum('ii',R)-1.0))
        skew = R - np.einsum('ij->ji',R)
        skew = np.array([skew[1,2],skew[0,2],-skew[0,1]]) 
        if abs(theta) > 1e-7:
            sin_term = 0.5 + (1.0/12.0) * theta ** 2 + (7.0/720.0) * theta ** 4 * (31.0/30240.0) * theta ** 6
        elif abs(theta) > 0.0:
            sin_term = theta / (2.0 * np.sin(theta))
        else:
            sin_term = 0.5
        return sin_term * skew
    
    def diff(self,X):
        temp_zero = np.zeros((3,1))
        return np.concatenate((X,temp_zero),axis=1) - np.concatenate((temp_zero,X),axis=1)

    def quad(self,X):
        temp_zero = np.zeros((3,1))
        return 0.5 * (np.concatenate((X,temp_zero),axis=1) + np.concatenate((temp_zero,X),axis=1))

    def update_kappa(self):
        temp_R = np.einsum('ijk,ilk->jlk',self.Q[:,:,1:], self.Q[:,:,:self.N-1])
        for ii in range(self.N-1):
            self.kappa[:,ii] = self.logm(temp_R[:,:,ii]) / self.D_0[ii]

    def update_sigma(self):
        self.sigma = np.einsum('ijk,ik->jk',self.Q, self.l / self.l_0 - self.Q[2,:,:])  

    def update_v(self, dvdt, dt):
        self.v = self.v + dvdt * dt
        
    def update_w(self, dwdt, dt):
        self.w = self.w + dwdt * dt

    def update_x(self, dt):
        self.x = self.x + self.v * dt
        self.l = self.x[:,1:] - self.x[:,:self.N]
        self.l_norm = np.linalg.norm(self.l,axis=0)
        self.D = 0.5 * (self.l_norm[1:] + self.l_norm[:self.N-1])
        self.e = self.l_norm / self.l_0
        self.ee = self.D / self.D_0
        
    def update_Q(self, dt):
        temp_expm = np.zeros((3,3,self.N))
        for ii in range(self.N):
            temp_expm[:,:,ii] = self.expm(self.w[:,ii] * dt)
        self.Q = np.einsum('ijk,jlk->ilk',temp_expm,self.Q)

    def update_acceleration(self):
        self.update_sigma()
        self.update_kappa()

        n_temp = np.einsum('ij,j...->j...', self.S, self.sigma)
        Qn_temp = np.einsum('ijk,jk->ik',self.Q,n_temp) / self.e
        dvdt = self.diff(Qn_temp) / self.m

        tau_temp = np.einsum('ij,j...->i...', self.B, self.kappa) / np.power(self.ee,3)
        kappa_temp = np.cross(self.kappa, tau_temp, axisa=0, axisb=0, axisc=0) * self.D_0
        shear_temp = np.einsum('ijk,jk->ik', self.Q, self.l / self.l_norm)
        shear_temp = np.cross(shear_temp, n_temp, axisa=0, axisb=0, axisc=0) * self.l_0
        dilatation_temp = np.einsum('ij,j...->i...',self.I,self.w) / self.e
        rigid_temp = np.cross(dilatation_temp, self.w, axisa=0, axisb=0, axisc=0)
        dilatation_temp /= self.e
        dwdt = self.diff(tau_temp) + self.quad(kappa_temp) + shear_temp + dilatation_temp + rigid_temp
        dwdt *= self.e
        dwdt = np.einsum('ij->ji',np.einsum('ij->ji', dwdt) * np.diag(self.I))

        return dvdt, dwdt

    def symplectic(self,timespan=10,dt=0.01,method='PEFRL'):
        if method == 'PEFRL':
            xi = 0.1786178958448091
            lmbda = -0.2123418310626054
            chi = -0.06626458266981849
            a = np.array([xi,chi,1-2*(xi+chi),chi,xi])
            b = np.array([0.5*(1-2*lmbda),lmbda,lmbda,0.5*(1-2*lmbda)])
        elif method == 'VV':
            a = np.array([0.5,0.5])
            b = np.array([1])
        elif method == 'PV':
            a = np.array([0,1])
            b = np.array([0.5,0.5])
        else:
            raise ValueError('Incompatible integrator specified')

        numsteps = int(timespan / dt)
        for ii in range(numsteps):
            for jj in range(max(len(a),len(b))):
                dvdt, dwdt = self.update_acceleration()
                if a[jj] != 0:
                    self.update_v(dvdt, a[jj] * dt)
                    self.update_w(dwdt, a[jj] * dt)
                if jj < len(b):
                    self.update_x(b[jj] * dt)
                    self.update_Q(b[jj] * dt)
            


filament = CosseratRod()
# print(filament.x)
# print(filament.Q)
# print(filament.B)
# print(filament.x[:,3])
# print(filament.expm(filament.x[:,0]))
# filament.update_Q(0.001)
# print(filament.logm(filament.Q[:,:,3]))
# print(filament.sigma)
# filament.update_sigma()
# print(filament.sigma)
# filament.update_kappa()
# print(filament.update_acceleration())
filament.symplectic()
print(filament.x)



