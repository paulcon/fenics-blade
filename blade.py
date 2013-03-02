from dolfin import *

def boundary(x,on_boundary):
    return on_boundary

if __name__=="__main__":
    mesh = Mesh("naca0018.xml")
    V = FunctionSpace(mesh,"CG",1)
    bc = DirichletBC(V,Constant(0.0),boundary)

    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(1.0)
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    plot(u,interactive=True)
