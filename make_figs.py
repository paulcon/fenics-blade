"""
Just a test file.
"""
from dolfin import *
from blade_solvers import *
import numpy as np

def load_kl_bases(V,d):
    bases = []
    for i in range(d):
        b = Function(V)
        File("kl/bases/%0.4d.xml.gz" % i) >> b
        bases.append(b)
        if not np.mod(i,50):
            info("Loaded basis %d" % i)
    return bases

if __name__ == "__main__":

    parameters["std_out_all_processes"] = False
    d = 250
    mesh = Mesh("mesh/naca0018_3d.xml.gz")
    V = FunctionSpace(mesh,"CG",1)
    bs = BladeSolver(V,d)
 
    a = np.random.randn(d)
    #a = np.zeros(d)
    bs.compute_q(a,True,'fig0') 
 
    bases = load_kl_bases(V,d)
    bs.compute_gradq(bases,True,'fig0')
   
