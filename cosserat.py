import numpy as np

class CosseratRod:

    def __init__(self, segments=10, seg_length=10.0, position=[[0.0],[0.0],[0.0]]):
        #Basic geometric properties
        self.N = segments
        self.l_0 = seg_length
        self.A = 38.5
        self.rho = 0.01

        #Current coordinates, velocities, strains, Voronoi, extensions,accelerations
        self.x = np.zeros((3,self.N+1))
        self.x[2,:] = np.array([i*self.l_0 for i in range(self.N+1)])
        self.x += np.array(position)
        self.v = np.zeros((3,self.N+1))
        self.Q = np.zeros((3,3,self.N))
        for ii in range(self.N):
            self.Q[:,:,ii] = np.array([[1,0,0],[0,1,0],[0,0,1]])
        self.w = np.zeros((3,self.N))
        self.kappa = np.zeros((3,self.N-1))
        self.sigma = np.zeros((3,self.N))
        self.l = self.x[:,1:] - self.x[:,:self.N]
        self.l_norm = np.linalg.norm(self.l,axis=0)
        self.D_0 = 0.5 * (self.l_norm[1:] + self.l_norm[:self.N-1])
        self.D = np.copy(self.D_0)
        self.e = self.l_norm / self.l_0
        self.ee = self.D / self.D_0
        self.e_v = np.einsum('ij,ij->j',self.l,self.v[:,1:] - self.v[:,:self.N]) / (self.l_norm * self.l_0);
        self.dvdt = np.zeros((3,self.N+1))
        self.dwdt = np.zeros((3,self.N))

        #Mass, inertia, rigidities
        self.m = (seg_length * self.A * self.rho) * np.ones((1,self.N+1))
        self.m[0,0] *= 0.5
        self.m[0,self.N] *= 0.5
        self.I = ((self.m[0,1] * self.A ** 2) / (4 * np.pi)) * np.diag(np.array([1.0,1.0,2.0]))
        self.B = np.diag(np.array([39000.0, 39000.0, 13000.0])) / 4.114
        self.S = np.diag(np.array([16000.0, 16000.0, 34000.0])) / 4.114
        
        #external forces
        self.ext_F = np.zeros((3,self.N+1))
        self.ext_C = np.zeros((3,self.N))

    def expm(self,w,theta):
        u = np.copy(w)
        #find norm of velocity vector u
        u_norm = np.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
        #check to see if velocity is close to zero
        if u_norm > 1e-14:
            #normalise velocity
            u /= u_norm
            #apply ax operator to convert vector into antisym matrix
            U = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
            #scale rotation angle by scalar velocity
            theta *= u_norm
            #implement Rodrigues formula
            return np.array([[1,0,0],[0,1,0],[0,0,1]]) + np.sin(theta) * U + (1 - np.cos(theta)) * np.einsum('ij,jk->ik',U, U)
        else:
            #if velocity is zero, no change in orientation
            return np.array([[1,0,0],[0,1,0],[0,0,1]])
    
    def logm(self,RR):
        R = np.copy(RR)
        #take arccos of interval (trace) to get angle theta
        theta = np.arccos(0.5 * (np.einsum('ii',R)-1.0) - 1e-10)
        #find different between R and R^T
        skew = R - np.einsum('ij->ji',R)
        #apply skew operator to convert antisym matrix to vector
        skew = np.array([skew[2,1],-skew[2,0],skew[1,0]])
        #conditionals to depending on value of angle
        if theta == 0:
            #no difference in orientation implies no curvature
            sin_term = 0
        elif abs(theta) > 1e-10:
            #analytic term for larg enough values
            sin_term = 0.5 * theta / (np.sin(theta))
        else:
            #taylor expansion for small values
            sin_term = 0.5 + (1.0/12.0) * (theta ** 2) + (7.0/720.0) * (theta ** 4) * (31.0/30240.0) * (theta ** 6)
        #return final vector term
        return -sin_term * skew
    
    def diff(self,X):
        #preallocate (3,1) vector of zeros 
        temp_zero = np.zeros((3,1))
        #implement difference operator
        return np.concatenate((X,temp_zero),axis=1) - np.concatenate((temp_zero,X),axis=1)

    def quad(self,X):
        #preallocate (3,1) vector of zeros
        temp_zero = np.zeros((3,1))
        #implement quadrature operator
        return 0.5 * (np.concatenate((X,temp_zero),axis=1) + np.concatenate((temp_zero,X),axis=1))

    def update_kappa(self):
        #calculate Q_i+1 * Q_i^T
        temp_R = np.einsum('ijk,ljk->ilk',self.Q[:,:,1:], self.Q[:,:,:self.N-1])
        #loop over kappas
        for ii in range(self.N-1):
            #calculate curvature between each segment
            self.kappa[:,ii] = self.logm(temp_R[:,:,ii]) / self.D_0[ii]

    def update_sigma(self):
        #implement sigma using einsum, equivalent to looped version
        self.sigma = np.einsum('ijk,jk->ik', self.Q, self.l / self.l_0 - self.Q[2,:,:])

    def update_v(self, dt, dissipation=0):
        #update filament linear velocity using dvdt, under dissipation
        self.v = self.v * (1 - dissipation * dt) + self.dvdt * dt
        #update velocity of dilatation
        self.e_v = np.einsum('ij,ij->j',self.l,self.v[:,1:] - self.v[:,:self.N]) / (self.l_norm * self.l_0);
        
    def update_w(self, dt, dissipation=0):
        #update filament angular velocity using dvdt, under dissipation
        self.w = self.w * (1 - dissipation * dt) + self.dwdt * dt

    def update_x(self, dt):
        #update position
        self.x = self.x + self.v * dt
        #update segment length vector
        self.l = self.x[:,1:] - self.x[:,:self.N]
        #calculate current segment length
        self.l_norm = np.sqrt(np.einsum('ij,ij->j',self.l, self.l))
        #calculate current Voronoi region length
        self.D = 0.5 * (self.l_norm[1:] + self.l_norm[:self.N-1])
        #calculate segment extension
        self.e = self.l_norm / self.l_0
        #calculate Voronoi extension
        self.ee = self.D / self.D_0
        
    def update_Q(self, dt):
        #preallocate rotation matrices
        temp_expm = np.zeros((3,3,self.N))
        #loop over orientation bases
        for ii in range(self.N):
            #calculate rotation matrices using expm
            temp_expm[:,:,ii] = self.expm(-self.w[:,ii], dt)
        #update orientations using rotation matrices
        self.Q = np.einsum('ijk,jlk->ilk',temp_expm,self.Q)

    def update_acceleration(self):
        #update sigma and kappa
        self.update_sigma()
        self.update_kappa()
        #calculate shear/stretch strain
        n_temp = np.einsum('ij,j...->j...', self.S, self.sigma)
        #take strain to laboratory frame
        Qn_temp = np.einsum('jik,jk->ik',self.Q,n_temp) / self.e 
        #update dvdt
        self.dvdt = (self.diff(Qn_temp) + self.ext_F) / self.m
        #calculate bend/twist strain
        tau_temp = np.einsum('ij,j...->i...', self.B, self.kappa) / np.power(self.ee,3)
        kappa_temp = np.cross(self.kappa, tau_temp, axisa=0, axisb=0, axisc=0) * self.D_0
        #calculate torsional effect of shear
        shear_temp = np.einsum('ijk,jk->ik', self.Q, self.l / self.l_norm)
        shear_temp = np.cross(shear_temp, n_temp, axisa=0, axisb=0, axisc=0) * self.l_0
        #temporarily define quantity
        dilatation_temp = np.einsum('ij,j...->i...',self.I,self.w) / self.e
        #calculate rigid body term (Euler)
        rigid_temp = np.cross(dilatation_temp, self.w, axisa=0, axisb=0, axisc=0)
        #calculate dilatation term
        dilatation_temp *= self.e_v / self.e
        #update dwdt
        self.dwdt = self.diff(tau_temp) + self.quad(kappa_temp) + shear_temp + dilatation_temp + rigid_temp + self.ext_C
        self.dwdt *= np.tile(self.e,(3,1)) / np.einsum('ij->ji',np.tile(np.diag(self.I),(10,1)))

    def update_conditions(self,conditions,numsteps,ii,dt):
        #no conditions set
        if conditions == []:
            return
        #clamped at 0 end
        if 'clamp_0' in conditions:
            self.v[:,:2] = np.zeros((3,2))
            self.w[:,0] = np.zeros((3,))
        #clamped at N end
        if 'clamp_N' in conditions:
            self.v[:,self.N-1:self.N+1] = np.zeros((3,2))
            self.w[:,self.N-1] = np.zeros((3,))
        #compressed and clamped at 0 end
        if 'comp_clamp_0' in conditions:
            self.v[:,:2] = np.zeros((3,2))
            if ii < numsteps // 10:
                self.v[2,:2] += 2.0 * 10 / (numsteps * dt)
        #compressed and clamped at N end
        if 'comp_clamp_N' in conditions:
            self.v[:,self.N-1:self.N+1] = np.zeros((3,2))
            if ii < numsteps // 10:
                self.v[2,self.N-1:self.N+1] += -2.0 * 10 / (numsteps * dt)
        #twisted and clamped at 0 end
        if 'twist_clamp_0' in conditions:
            self.w[:,0] = np.zeros((3,))
            if ii < numsteps // 10:
                self.w[2,0] += 1.0 * 10 / (numsteps * dt)
        #twisted and clamped at N end
        if 'twist_clamp_N' in conditions:
            self.w[:,self.N-1] = np.zeros((3,))
            if ii < numsteps // 10:
                self.w[2,self.N-1] += -1.0 * 10 / (numsteps * dt)
        #add noise to simulations
        if 'noise_no_ext' in conditions:
            force_scale = 0.01
            self.ext_F = np.random.normal(size=(3,self.N+1)) * force_scale
            self.ext_C = np.random.normal(size=(3,self.N)) * force_scale



    def symplectic(self,timespan=100,dt=0.01,method='PEFRL',ext_F=None,ext_C=None,dissipation=0,conditions=[]):
        #Check to see if there are any external forces acting on filament, else use zero arrays
        if isinstance(ext_F,np.ndarray) == False:
            self.ext_F = np.zeros((3,self.N+1))
        else:
            self.ext_F = ext_F
        if isinstance(ext_C,np.ndarray) == False:
            self.ext_C = np.zeros((3,self.N))
        else:
            self.ext_C = ext_C

        #Set integrating method as position extended Forrest-Ruth like   
        if method == 'PEFRL':
            xi = 0.1786178958448091
            lmbda = -0.2123418310626054
            chi = -0.06626458266981849
            a = np.array([xi,chi,1-2*(xi+chi),chi,xi])
            b = np.array([0.5*(1-2*lmbda),lmbda,lmbda,0.5*(1-2*lmbda)])
        #Set integrating method as velocity Verlet
        elif method == 'VV':
            a = np.array([0.5,0.5])
            b = np.array([1])
        #Set integrating method as position Verlet
        elif method == 'PV':
            a = np.array([0,1])
            b = np.array([0.5,0.5])
        #Raise error if method not recognised
        else:
            raise ValueError('Incompatible integrator specified')

        #Set number of steps in simulation
        numsteps = int(timespan / dt)
        #Run simulation loop
        for ii in range(numsteps):
            #Loop through integrator coefficents
            for jj in range(max(len(a),len(b))):
                #check to see if there is a first velocity step
                if a[jj] != 0:
                    #update accelerations
                    self.update_acceleration()
                    #update velocities
                    self.update_v(a[jj] * dt, dissipation)
                    self.update_w(a[jj] * dt, dissipation)
                    #enforce boundary conditions
                    self.update_conditions(conditions,numsteps,ii,dt)
                #check to see if there is a final position step
                if jj < len(b):
                    #update positions and orientations
                    self.update_x(b[jj] * dt)
                    self.update_Q(b[jj] * dt)


