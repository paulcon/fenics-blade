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

   d = 10
   mesh = Mesh("mesh/naca0018_3d.xml.gz")
   V = FunctionSpace(mesh,"CG",1)
   bs = BladeSolver(V,d)
   
   a = np.zeros(d)
   q0 = bs.compute_q(a)
   klbases = load_kl_bases(V,d)
   da_adj = bs.compute_gradq(klbases) 
   
   eps = 1e-6
   da_fd = np.zeros(d)
   for i in range(d):
      aa = np.zeros(d); aa[i]=eps
      q1 = bs.compute_q(aa)
      da_fd[i] = (1/eps)*(q1-q0)
   pdb.set_trace()
      
      
   
