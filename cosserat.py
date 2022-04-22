import numpy as np

class CosseratRod:

    def __init__(self, segments=10, seg_length=10, position=[[0],[0],[0]]):
        #Basic geometric properties
        self.N = segments
        self.l_0 = seg_length
        self.A = 38.5
        self.rho = 0.01

        #Current coordinates, velocities, strains
        self.x = np.zeros((3,self.N+1))
        self.x[0,:] = np.array([i*self.l_0 for i in range(self.N+1)])
        self.x += np.array(position)
        self.v = np.zeros((3,self.N+1))
        self.Q = np.tile(np.identity(3), self.N)
        self.w = np.zeros((3,self.N))
        self.kappa = np.zeros((3,self.N-1))
        self.sigma = np.zeros((3,self.N))

        #Mass, inertia, rigidities
        self.m = (seg_length * self.A * self.rho) * np.ones((1,self.N+1))
        self.m[0,0] /= 2
        self.m[0,self.N] /= 2
        self.I = ((self.m[0,1] * self.A ** 2) / (4 * np.pi)) * np.diag(np.array([1,1,2]))
        self.B = np.diag(np.array([39000, 39000, 13000]))
        self.S = np.diag(np.array([16000, 16000, 34000]))
        
    def expm(self,u):
        theta = np.linalg.norm(u)
        if theta > 0:
            u /= theta
            U = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
            return np.identity(3) + np.sin(theta) * U + np.cos(theta) * U * U
        return np.identity(3)
    
    def logm(self,R):
        pass
    
    def update_kappa(self):
        pass

    def update_sigma(self):
        pass

    def update_v(self, dvdt, dt):
        self.v = self.v + dvdt * dt
        
    def update_w(self, dwdt, dt):
        self.w = self.w + dwdt * dt

    def update_x(self, dt):
        self.x = self.x + self.v * dt
        
    def update_Q(self, dt):
        for ii in range(self.N):
            self.Q[:,3*ii:3*(ii+1)] = np.einsum('ij,jk->ik', self.expm(self.w[:,ii] * dt), self.Q[:,3*ii:3*(ii+1)])

    def update_acceleration(self):
        pass

    def symplectic(self,timespan,dt,coeffs):
        pass


filament = CosseratRod()
print(filament.x)
print(filament.Q)
print(filament.B)
print(filament.x[:,3])
print(filament.expm(filament.x[:,0]))
filament.update_Q(0.001)


