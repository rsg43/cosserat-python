import numpy as np

class CosseratRod:

    def __init__(self, segments=10, seg_length=10.0, position=[[0.0],[0.0],[0.0]]):
        #Basic geometric properties
        self.N = segments
        self.l_0 = seg_length
        self.A = 38.5
        self.rho = 0.01

        #Current coordinates, velocities, strains
        self.x = np.zeros((3,self.N+1))
        self.x[2,:] = np.array([i*self.l_0 for i in range(self.N+1)])
        self.x += np.array(position)
        self.v = np.zeros((3,self.N+1))
        # self.Q = np.tile(np.identity(3), self.N)
        self.Q = np.zeros((3,3,self.N))
        for ii in range(self.N):
            self.Q[:,:,ii] = np.identity(3)
        self.w = np.zeros((3,self.N))
        self.kappa = np.zeros((3,self.N-1))
        self.sigma = np.zeros((3,self.N))
        self.l = self.x[:,1:] - self.x[:,:self.N]
        print(self.l)
 

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
        #Change this to act on a 3xN array of u values
    
    def logm(self,R):
        theta = np.arccos(0.5 * (np.einsum('ii',R)-1.0))
        skew = R - np.einsum('ij->ji',R)
        skew = np.array([[skew[1,2]],[skew[0,2]],[-skew[0,1]]]) 
        if abs(theta) > 1e-7:
            sin_term = 0.5 + (1.0/12.0) * theta ** 2 + (7.0/720.0) * theta ** 4 * (31.0/30240.0) * theta ** 6
        elif abs(theta) > 0.0:
            sin_term = theta / (2.0 * np.sin(theta))
        else:
            sin_term = 0.5
        return sin_term * skew
        #write this to act on 3x3N array of orientation bases
    
    def update_kappa(self):
        pass

    def update_sigma(self):
        self.sigma = np.einsum('ijk...,ik...->jk...',self.Q, self.l / self.l_0 - self.Q[2,:,:])  

    def update_v(self, dvdt, dt):
        self.v = self.v + dvdt * dt
        
    def update_w(self, dwdt, dt):
        self.w = self.w + dwdt * dt

    def update_x(self, dt):
        self.x = self.x + self.v * dt
        
    def update_Q(self, dt):
        temp_expm = np.zeros((3,3,self.N))
        for ii in range(self.N):
            temp_expm[:,:,ii] = self.expm(self.w[:,ii] * dt)
        #try to get rid of this for loop if possible
        self.Q = np.einsum('ijk,jlk->ilk',temp_expm,self.Q)

    def update_acceleration(self):
        return dvdt, dwdt

    def symplectic(self,timespan,dt,coeffs):
        pass


filament = CosseratRod()
print(filament.x)
print(filament.Q)
print(filament.B)
print(filament.x[:,3])
print(filament.expm(filament.x[:,0]))
filament.update_Q(0.001)
print(filament.logm(filament.Q[:,:,3]))
# print(filament.sigma)
# filament.update_sigma()
# print(filament.sigma)


