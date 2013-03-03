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
        self.bases = np.loadtxt(open("klbases_2d.txt","rb"),delimiter=",")
        self.bmeshx = np.loadtxt(open("xcoord.txt","rb"),delimiter=",")
        self.bmeshz = np.loadtxt(open("zcoord.txt","rb"),delimiter=",")

    def construct(self,coeff):
        # linear combination of bases and parameters
        rf = np.dot(self.bases,coeff)
        self.f = RectBivariateSpline(self.bmeshx,self.bmeshz, \
                                     rf.reshape(len(self.bmeshx), \
                                                len(self.bmeshz)))

    def eval(self,values,x):
        values[0] = self.f(x[0],x[2])

if __name__=="__main__":

    # Load the mesh and define the function space
    mesh = Mesh("naca0018_3d.xml")
    V = FunctionSpace(mesh,"CG",1)

    # Dirichlet boundary
    bc_cool = CoolBoundary()
    bcs = DirichletBC(V,Constant(0.0),bc_cool)

    # Random Neumann boundary
    f = RandomHeatFlux()
    a = np.random.normal(0,1.0,500);
    f.construct(a)

    # Write heat flux
    fluxfile = File("flux3d.pvd")
    fluxfile << project(f,V)
    
    # Set up
    u = TrialFunction(V)
    v = TestFunction(V)
    kappa = Constant(1.0) + Constant(0.01)*u
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

    solnfile = File("blade3d.pvd")
    solnfile << u_
    
    #plot(u_,interactive=True)
