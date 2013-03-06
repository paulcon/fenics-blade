"""
Load the text files for the KL bases, and construct DOLFIN objects.
Then save them to files.
"""
from dolfin import *
from scipy.interpolate import RectBivariateSpline
import numpy as np

class KLBasis(Expression):
    def __init__(self,bmeshx,bmeshz,b):
        self.f = RectBivariateSpline(bmeshx,bmeshz, \
                                     b.reshape(len(bmeshx), len(bmeshz)))
    def eval(self,values,x):
        tol = 1E-14 # tolerance for coordinate comparisons
        if abs(x[2])<tol or abs(x[2]-2)<tol:
            values[0] = 0.0
        else:
            values[0] = self.f(x[0],x[2])

if __name__ == "__main__":

    # load mesh and build function space
    mesh = Mesh("../mesh/naca0018_3d.xml.gz")
    V = FunctionSpace(mesh,"CG",1)
    
    bases = np.loadtxt(open("klbases_2d.txt","rb"),delimiter=",")
    bmeshx = np.loadtxt(open("xcoord.txt","rb"),delimiter=",")
    bmeshz = np.loadtxt(open("zcoord.txt","rb"),delimiter=",")

    for i in range(500):
        b = bases[:,i]
        klb = KLBasis(bmeshx,bmeshz,b)
        bfilename = "bases/%0.4d.xml.gz" % i
        bfile = File(bfilename)
        bfile << project(klb,V)
