from dolfin import *
import numpy as np

class CoolBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        circ1 = (x[0]-0.2)*(x[0]-0.2) + x[1]*x[1] - 0.0025 < tol
        circ2 = (x[0]-0.45)*(x[0]-0.45) + x[1]*x[1] - 0.0025 < tol
        circ3 = (x[0]-0.7)*(x[0]-0.7) + x[1]*x[1] - 0.0009 < tol
        return on_boundary and (circ1 or circ2 or circ3)
    

class HotFluxBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        top = 0.5*(0.2969*np.sqrt(x[0]) - 0.1260*x[0] \
                   - 0.3516*x[0]*x[0] + 0.2843*x[0]*x[0]*x[0] \
                   - 0.1015*x[0]*x[0]*x[0]*x[0]) - x[1] < tol
        bot = 0.5*(0.2969*np.sqrt(x[0]) - 0.1260*x[0] \
                   - 0.3516*x[0]*x[0] + 0.2843*x[0]*x[0]*x[0] \
                   - 0.1015*x[0]*x[0]*x[0]*x[0]) + x[1] < tol
        return on_boundary and (top or bot)

if __name__=="__main__":
    mesh = Mesh("naca0018_2d.xml")
    V = FunctionSpace(mesh,"CG",2)

    #boundary_parts = \
    #    MeshFunction("uint", mesh, mesh.topology().dim()-1)

    
    bc_cool = CoolBoundary()
    bcs = DirichletBC(V,Constant(0.0),bc_cool)
    
    
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(1.0)
    g = Expression("1-x[0]")
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx + g*v*ds

    # Compute solution
    u = Function(V)
    solve(a == L, u, bcs)

    plot(u,interactive=True)
