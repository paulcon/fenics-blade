"""
This was an attempt to use dolfin-adjoint, but it didn't work very well.

There is a way to fix this. It requires building the heat flux function
from DOLFIN objects. But it seems a bit more trouble than it's worth
for now -- especially since we can hand code the adjoint solver just
as easily as the forward solver.
"""
from dolfin import *
from dolfin_adjoint import *
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
        rf = np.dot(self.bases,coeff)
        self.f = RectBivariateSpline(self.bmeshx,self.bmeshz, \
                                     rf.reshape(len(self.bmeshx), \
                                                len(self.bmeshz)))

    def eval(self,values,x):
        values[0] = self.f(x[0],x[2])

if __name__=="__main__":
    #pdb.set_trace()
    
    mesh = Mesh("mesh/naca0018_3d.xml")
    V = FunctionSpace(mesh,"CG",1)

    # Dirichlet boundary
    bc_cool = CoolBoundary()
    bcs = DirichletBC(V,Constant(0.0),bc_cool)

    # Random Neumann boundary
    f = RandomHeatFlux()
    a = np.random.normal(0,1.0,500);
    f.construct(a)
    fV = project(f,V)

    # Write heat flux
    fluxfile = File("output/flux3d.pvd")
    fluxfile << fV
    
    # Set up
    u = TrialFunction(V)
    v = TestFunction(V)
    kappa = Constant(1.0) + Constant(0.01)*u
    F = -inner(kappa*grad(u), grad(v))*dx + fV*v*ds

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

    pdb.set_trace()
    g = Expression("x[0]>0.85")
    q = Functional(inner(g,u_)*dx)

    # Yeah, this doesn't work.
    dqda = compute_gradient(q,InitialConditionParameter(fV))

    solnfile = File("output/blade3d.pvd")
    solnfile << u_
    
    #plot(u_,interactive=True)
