import mk
from scipy.sparse.linalg import bicg

class solver:
    @staticmethod
    def solve_dU(K,Udof,dUb,Fdof,dFb):
        ERR = 10**(-6)
        Kb, rhs = mk.mk.mk_global_mtr(K, Udof, dUb, Fdof, dFb)
        dU,FLG = bicg(Kb,rhs,tol=ERR,maxiter=5000)
        return dU.reshape(-1,1),FLG
