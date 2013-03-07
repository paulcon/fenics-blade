"""
Just a test file.
"""
from dolfin import *
from blade_solvers import *
import numpy as np

if __name__ == "__main__":

   d = 10
   mesh = Mesh("mesh/naca0018_3d.xml.gz")
   V = FunctionSpace(mesh,"CG",1)
   bs = BladeSolver(V,d)

   a = np.random.rand(d)
   bs.compute_q(a,True,'test')
   
