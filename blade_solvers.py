from dolfin import *
import numpy as np
from scipy.interpolate import RectBivariateSpline
import pdb


class CoolBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        circ1 = (x[0]-0.2)*(x[0]-0.2) + x[1]*x[1] - 0.0025 < tol
        circ2 = (x[0]-0.45)*(x[0]-0.45) + x[1]*x[1] - 0.0025 < tol
        circ3 = (x[0]-0.7)*(x[0]-0.7) + x[1]*x[1] - 0.0009 < tol
        return on_boundary and (circ1 or circ2 or circ3)

class RandomHeatFlux(Expression):
    def __init__(self):
        # load 2d bases and mesh they live on
        self.bases = np.loadtxt(open("kl/klbases_2d.txt","rb"),delimiter=",")
        self.bmeshx = np.loadtxt(open("kl/xcoord.txt","rb"),delimiter=",")
        self.bmeshz = np.loadtxt(open("kl/zcoord.txt","rb"),delimiter=",")
        
    def construct(self,coeff):
        # linear combination of bases and parameters
        d = len(coeff)
        rf = np.dot(self.bases[:,:d],coeff)
        self.f = RectBivariateSpline(self.bmeshx,self.bmeshz, \
                                     rf.reshape(len(self.bmeshx), \
                                                len(self.bmeshz)))

    def eval(self,values,x):
        tol = 1E-14 # tolerance for coordinate comparisons
        if abs(x[2])<tol or abs(x[2]-2)<tol:
            values[0] = 0.0
        else:
            # model for random heat flux
            mean = (max(650.0-1950.0*x[0],0.0) + 250.0*x[0]*(np.tanh(7.0*x[0]-4.0)+1.0) + 300.0)/1000.0
            sig = (100.0*np.exp(-10.0*(x[0]-0.8)**2)+5)/1000.0
            values[0] = mean + sig*self.f(x[0],x[2])

class BladeSolver():
    def __init__(self,V,d):

        self.V = V 
        self.d = d 

        # conductivity parameters
        self.k0 = 1.0
        self.k1 = 0.01
        
        # boundary conditions
        bc_cool_fwd = CoolBoundary()
        self.bcs_fwd = DirichletBC(self.V,Constant(0.35),bc_cool_fwd)
        bc_cool_adj = CoolBoundary()
        self.bcs_adj = DirichletBC(self.V,Constant(0.0),bc_cool_adj)

        # construct heat flux
        self.heat_flux = RandomHeatFlux()

        # inner product for QoI and rhs for adjoint solve
        self.g = Expression("x[0]>0.3")
        
    def compute_q(self,coeff,savesoln=False,filename=None):
        self.heat_flux.construct(coeff)

        # Set up
        u = TrialFunction(self.V)
        v = TestFunction(self.V)
        kappa0 = Constant(self.k0)
        kappa1 = Constant(self.k1)
        kappa = kappa0 + kappa1*u
        F = -inner(kappa*grad(u), grad(v))*dx + self.heat_flux*v*ds

        # Taking derivative
        u_ = Function(self.V)
        F = action(F,u_)
        J = derivative(F,u_,u)

        # Nonlinear solver
        fwd_problem = NonlinearVariationalProblem(F, u_, self.bcs_fwd, J)
        fwd_solver = NonlinearVariationalSolver(fwd_problem)
        prm = fwd_solver.parameters
        prm["nonlinear_solver"] = "snes"
        prm["snes_solver"]["maximum_iterations"] = 25
        prm["snes_solver"]["report"] = False
        #pdb.set_trace()
        fwd_solver.solve()

        if savesoln:
            solnfile = File("output/temp_%s_.pvd" % filename)
            solnfile << u_
            fluxfile = File("output/flux_%s_.pvd" % filename)
            fluxfile << project(self.heat_flux,self.V)

        self.u_ = u_
        return assemble(inner(self.g,u_)*dx)

    def compute_gradq(self,klbases,savesoln=False,filename=None):
        
        w = TrialFunction(self.V)
        v = TestFunction(self.V)
        kappa0 = Constant(self.k0)
        kappa1 = Constant(self.k1)
        a = inner(grad(v),(kappa0+kappa1*self.u_)*grad(w))*dx \
            - kappa1*v*dot(grad(self.u_),grad(w))*dx
        L = self.g*v*dx

        w_ = Function(self.V)
        adj_problem = LinearVariationalProblem(a,L,w_,self.bcs_adj)
        adj_solver = LinearVariationalSolver(adj_problem)
        adj_solver.solve()

        if savesoln:
            adjfile = File("output/adj_%s_.pvd" % filename)
            adjfile << w_
        
        gradq = np.zeros(self.d)
        for i in range(self.d):
            gradq[i] = assemble(inner(w_,klbases[i])*ds)
        return gradq
        
