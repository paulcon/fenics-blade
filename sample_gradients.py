from dolfin import *
import numpy as np
from blade_solvers import *
import pdb
import time

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
    #pdb.set_trace()
    parameters["std_out_all_processes"] = False
    
    # load the design
    design =  np.loadtxt(open("design/initial_design.txt","rb"),delimiter=",")
    [n,d] = design.shape

    # construct the solver
    mesh = Mesh("mesh/naca0018_3d.xml.gz")
    V = FunctionSpace(mesh,"CG",1)
    bs = BladeSolver(V,d)

    # load KL basis as DOLFIN objects
    kl_bases = load_kl_bases(V,d)

    #pdb.set_trace()
    q = np.zeros(n)
    gradq = np.zeros([n,d])
    timez = np.zeros([n,2])
    for i in range(n):
        coeff = design[i,:]

        # compute QoI
        tt = time.clock()
        q[i] = bs.compute_q(coeff)
        timez[i,0] =  (time.clock() - tt)/60.0
        info("Time to compute q: %4.2f (min)" % timez[i,0])

        # compute adjoint and gradients
        tt = time.clock()
        gradq[i,:] = bs.compute_gradq(kl_bases)
        timez[i,1] =  (time.clock() - tt)/60.0
        info("Time to compute gradq: %4.2f (min)" % timez[i,1])
        
        info("QoI[%d] = %6.4f" % (i,q[i]))

    np.savetxt('output/qois.txt',q,delimiter=",",fmt="%18.16e")
    np.savetxt('output/gradients.txt',gradq,delimiter=",",fmt="%18.16e")
        
    

