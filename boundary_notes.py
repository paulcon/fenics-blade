tol = 1E-14
def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]) < tol
def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-1) < tol
Gamma_0 = DirichletBC(V, Constant(0.0), left_boundary)
Gamma_1 = DirichletBC(V, Constant(1.0), right_boundary)
bcs = [Gamma_0, Gamma_1]

# boudary parts
boundary_parts = \
MeshFunction("uint", mesh, mesh.topology().dim()-1)

class LowerRobinBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        return on_boundary and abs(x[1]) < tol
    
Gamma_R = LowerRobinBoundary()
Gamma_R.mark(boundary_parts, 0)

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        return on_boundary and abs(x[0]) < tol

Gamma_0 = LeftBoundary()
Gamma_0.mark(boundary_parts, 2)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14 # tolerance for coordinate comparisons
        return on_boundary and abs(x[0] - 1) < tol

Gamma_1 = RightBoundary()
Gamma_1.mark(boundary_parts, 3)
