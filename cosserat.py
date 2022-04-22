import numpy as np

class CosseratRod:

    def __init__(self, segments=10, seg_length=10, mass=1):
        self.N = segments
        self.l_0 = seg_length

        self.x = np.zeros((3,N+1))
        self.v = np.zeros((3,N+1))
        self.Q = [np.identity(3) for i in range(N)]
        self.w = np.zeros((3,N))
        self.kappa = np.zeros((3,N-1))
        self.sigma = np.zeros((3,N))

        self.m = mass * np.ones((1,N+1))
        self.n[0] /= 2
        self.n[N] /= 2

        self.I = np.identity(np.array([]))
        self.B = np.identity(np.array([]))
        self.S = np.identity(np.array([]))
        
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


filament = CosseratRod
print(filament)
