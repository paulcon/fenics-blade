"""
Solve the heat transfer on the blade with a random heat flux. More notes coming soon.
Include the adjoint solution.
"""

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
    def __init__(self,mean,sig):
        # load 2d bases and mesh they live on
        self.bases = np.loadtxt(open("kl/klbases_2d.txt","rb"),delimiter=",")
        self.bmeshx = np.loadtxt(open("kl/xcoord.txt","rb"),delimiter=",")
        self.bmeshz = np.loadtxt(open("kl/zcoord.txt","rb"),delimiter=",")
        self.mean = mean
        self.sig = sig
        
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
            values[0] = self.mean + self.sig*self.f(x[0],x[2])

def load_kl_bases(V,d):
    bases = []
    for i in range(d):
        b = Function(V)
        File("kl/bases/%0.4d.xml.gz" % i) >> b
        bases.append(b)
        print "Loaded basis %d" % i
    return bases
    
if __name__=="__main__":

    # number of parameters
    d = 500
    
    # Load the mesh and define the function space
    mesh = Mesh("mesh/naca0018_3d.xml.gz")
    V = FunctionSpace(mesh,"CG",1)

    # Dirichlet boundary
    bc_cool = CoolBoundary()
    bcs = DirichletBC(V,Constant(0.0),bc_cool)

    # Random Neumann boundary
    mean = 1.0; sig = 0.2
    f = RandomHeatFlux(mean,sig)
    a = np.random.normal(0,1.0,d);
    f.construct(a)

    # Write heat flux
    fluxfile = File("output/flux3d.pvd")
    fluxfile << project(f,V)
    
    # Set up
    u = TrialFunction(V)
    v = TestFunction(V)
    kappa0 = Constant(1.0)
    kappa1 = Constant(0.01)
    kappa = kappa0 + kappa1*u
    F = -inner(kappa*grad(u), grad(v))*dx + f*v*ds

    # Taking derivative
    u_ = Function(V)
    F = action(F,u_)
    J = derivative(F,u_,u)

    # Solve the nonlinear problem
    problem = NonlinearVariationalProblem(F, u_, bcs, J)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters
    prm["nonlinear_solver"] = "snes"
    prm["snes_solver"]["maximum_iterations"] = 25
    #set_log_level(PROGRESS)
    solver.solve()

    # Compute the quantity of interest
    g = Expression("x[0]>0.85")
    q = inner(g,u_)*dx
    print "QoI: %6.4e" % assemble(q)

    #pdb.set_trace()
    # Solve the adjoint problem
    #bc_cool_adj = CoolBoundary()
    #bcs_adj = DirichletBC(V,Constant(0.0),bc_cool)
    w_adj = TrialFunction(V)
    a_adj = inner(grad(v),(kappa0+kappa1*u_)*grad(w_adj))*dx - kappa1*v*dot(grad(u_),grad(w_adj))*dx
    L_adj = g*v*dx

    w = Function(V)
    solve(a_adj == L_adj, w, bcs)

    
    #solnfile = File("output/blade3d.pvd")
    #solnfile << u_
    bases = load_kl_bases(V,d)
    #pdb.set_trace()
    da = np.zeros(d)
    for i in range(d):
        da[i] = assemble(inner(w,bases[i])*ds)

    #pdb.set_trace()
    print "done!"

    
        
    
