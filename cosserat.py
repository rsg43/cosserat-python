import numpy as np

class Filament:
    def __init__(self, segments=10, L=10, m=1):
        self.m = m
        self.N = segements
        self.L = L
        self.x = np.zeros((3,N+1))
        self.v = np.zeros((3,N+1))
        self.Q = [np.identity(3) for i in range(N)]
        self.w = np.zeros((3,N))
        self.kappa = np.zeros((3,N-1))
        self.sigma = np.zeros((3,N))
        
    def expm(u):
        pass
    
    def logm(R):
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
        for i in range(self.N):
            self.Q[i] = expm(self.w[i] * dt) * self.Q[i]

    def update_acceleration(self):
        pass

    def symplectic(self,timespan,dt,coeffs):
        pass

