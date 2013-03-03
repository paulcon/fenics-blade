from dolfin import *
import numpy as np
#from scipy.interpolate import interp1d
import pdb

class CoolBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        circ1 = (x[0]-0.2)*(x[0]-0.2) + x[1]*x[1] - 0.0025 < tol
        circ2 = (x[0]-0.45)*(x[0]-0.45) + x[1]*x[1] - 0.0025 < tol
        circ3 = (x[0]-0.7)*(x[0]-0.7) + x[1]*x[1] - 0.0009 < tol
        return on_boundary and (circ1 or circ2 or circ3)
    
class TopFluxBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        top = 0.5*(0.2969*np.sqrt(x[0]) - 0.1260*x[0] \
                   - 0.3516*x[0]*x[0] + 0.2843*x[0]*x[0]*x[0] \
                   - 0.1015*x[0]*x[0]*x[0]*x[0]) - x[1] < tol
        return on_boundary and top

class BottomFluxBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        bot = 0.5*(0.2969*np.sqrt(x[0]) - 0.1260*x[0] \
                   - 0.3516*x[0]*x[0] + 0.2843*x[0]*x[0]*x[0] \
                   - 0.1015*x[0]*x[0]*x[0]*x[0]) + x[1] < tol
        return on_boundary and bot

class RandomHeatFlux(Expression):
    def __init__(self):
        # load 2d bases and mesh they live on
        self.bases = np.loadtxt(open("klbases_1d.txt","rb"),delimiter=",")
        self.bases_mesh = np.loadtxt(open("xcoord.txt","rb"),delimiter=",")

    def construct(self,coeff):
        # linear combination of bases and parameters
        #self.f = interp1d(self.bases_mesh,self.bases*coeff)
        self.f = 1.0 + np.dot(self.bases,coeff)

    def eval(self,values,x):
        #values[0] = ... function of (x[0],x[2])
        values[0] = np.interp(x[0],self.bases_mesh,self.f)

if __name__=="__main__":
    #pdb.set_trace()
    
    mesh = Mesh("naca0018_2d.xml")
    V = FunctionSpace(mesh,"CG",4)

    # Dirichlet boundary
    bc_cool = CoolBoundary()
    bcs = DirichletBC(V,Constant(0.0),bc_cool)

    # Random Neumann boundary
    f = RandomHeatFlux()
    a = np.random.normal(0,1.0,31);
    f.construct(a)

    u = TrialFunction(V)
    v = TestFunction(V)
    #g = Expression("1-10*(sin(4*3.1415*x[0])-0.3*cos(3.1415*x[0]))")
    kappa = Constant(1.0) + Constant(0.01)*u
    F = -inner(kappa*grad(u), grad(v))*dx + f*v*ds

    # Taking derivative
    u_ = Function(V)
    F = action(F,u_)
    J = derivative(F,u_,u)

    problem = NonlinearVariationalProblem(F, u_, bcs, J)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters
    prm["newton_solver"]["absolute_tolerance"] = 1E-8
    prm["newton_solver"]["relative_tolerance"] = 1E-7
    prm["newton_solver"]["maximum_iterations"] = 25

    set_log_level(PROGRESS)
    solver.solve()

    qoi = np.power(assemble((u_**8)*dx),1.0/8.0)
    print "QoI: %6.4e" % qoi

    plot(u_,interactive=True)
