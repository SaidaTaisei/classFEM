import numpy as np
from scipy.sparse import csr_matrix
import sys

class mk:
    @staticmethod
    def mk_stif_mtr(node,elem,x,y,z,U,dU,E,H,nu,S,smax):
        num = node.shape[0]
        enum = elem.shape[0]
        DOF = 3*num
        poi = smax.shape[0]
        itg = (poi/enum)**(1/3)

        K = np.zeros([DOF,DOF])
        dS = np.zeros([poi,6])
        sm = np.zeros([poi,1])

        for ee in range(1,enum+1):
            idx = elem[ee-1,:]
            dof = np.concatenate([idx,idx+num,idx+2*num])
            gpx = np.array(range(1,int(itg)**3+1)) + (ee-1)*int(itg)**3
            gpx = gpx.reshape(gpx.shape[0],1)

            Ee = E[ee-1]
            He = H[ee-1]
            nue = nu[ee-1]
            enode = node[idx-1,:]
            Se = S[np.ravel(gpx)-1,:]
            semax = smax[np.ravel(gpx)-1]

            xe = x[idx - 1]
            ye = y[idx - 1]
            ze = z[idx - 1]
            Ue = U[np.concatenate([idx - 1, idx - 1 + num, idx - 1 + 2 * num])].reshape(-1,1)
            dUe = dU[np.concatenate([idx - 1, idx - 1 + num, idx - 1 + 2 * num])].reshape(-1,1)
            Ke,dSe,sme = mk.mk_elem_stif_mtr(enode, xe, ye, ze, Ue, dUe, Ee, He, nue, Se, semax)
            #print(Ke.shape)
            #print(K[dof-1,dof-1].shape)
            #print(dof.shape)
            dS[np.ravel(gpx)-1, :] = dSe.T
            K[np.ix_(dof-1, dof-1)] = K[np.ix_(dof-1, dof-1)] + Ke
            smax[np.ravel(gpx)-1] = np.maximum(sme,semax)
        sMax = np.copy(smax)
        return K,dS,sMax

    @staticmethod
    def mk_elem_stif_mtr(node,x,y,z,U,dU,E,H,nu,Sij,smax):
        poi = smax.shape[0]
        itg = int(poi ** (1 / 3))
        P, W = mk.lgwt(itg,-1,1)
        X = node[:, 3].reshape(-1,1)
        Y = node[:, 4].reshape(-1,1)
        Z = node[:, 5].reshape(-1,1)
        u = U[0:20, 0].reshape(-1,1)
        v = U[20:40, 0].reshape(-1,1)
        w = U[40:60, 0].reshape(-1,1)
        du = dU[0:20, 0].reshape(-1,1)
        dv = dU[20:40, 0].reshape(-1,1)
        dw = dU[40:60, 0].reshape(-1,1)
        u = u + du
        v = v + dv
        w = w + dw
        x = x + du
        y = y + dv
        z = z + dw
        lbd0 = nu * E / ((1 + nu) * (1 - 2 * nu))
        mu0 = E / (2 * (1 + nu))
        lbd1 = nu * H / ((1 + nu) * (1 - 2 * nu))
        mu1 = H / (2 * (1 + nu))
        D0 = mk.mk_D_mtr(lbd0, mu0)
        D1 = mk.mk_D_mtr(lbd1, mu1)
        cc = 0
        dSij = np.zeros([6, poi])
        ms = np.zeros([poi, 1])
        Ke = np.zeros([60, 60])
        for xx in range(1,itg+1):
            for yy in range(1, itg+1):
                for zz in range(1, itg+1):
                    cc = cc + 1
                    weight = W[xx-1] * W[yy-1] * W[zz-1]
                    B, F, J = mk.mk_B_mtr(P[xx-1], P[yy-1], P[zz-1], X, Y, Z, x, y, z, u, v, w)
                    dEij = np.dot(B.T, dU)
                    ms[cc-1] = mk.mk_Mises_stress(Sij[cc-1,:], F)
                    if(ms[cc-1]**2 >= smax[cc-1]**2):
                        dSij[:, cc-1] = np.dot(D1, dEij)[:,0]
                        Ke = Ke + weight * J * np.dot(np.dot(B, D1), B.T)
                    else:
                        dSij[:, cc-1] = np.dot(D0, dEij)[:,0]
                        Ke = Ke + weight * J * np.dot(np.dot(B, D0), B.T)
        return Ke, dSij, ms

    @staticmethod
    def mk_D_mtr(lbd,mu):
        D = np.zeros([6,6])
        D[0, 0] = lbd + 2 * mu
        D[0, 1] = lbd
        D[0, 2] = lbd
        D[1, 0] = lbd
        D[1, 1] = lbd + 2 * mu
        D[1, 2] = lbd
        D[2, 0] = lbd
        D[2, 1] = lbd
        D[2, 2] = lbd + 2 * mu
        D[3, 3] = mu
        D[4, 4] = mu
        D[5, 5] = mu
        return D

    @staticmethod
    def mk_B_mtr(P,Q,R,X,Y,Z,x,y,z,u,v,w):
        dNdP = mk.mk_dNdP_vec(P, Q, R, 1)
        dNdQ = mk.mk_dNdP_vec(P, Q, R, 2)
        dNdR = mk.mk_dNdP_vec(P, Q, R, 3)
        #print("mkB")
        #print(dNdP.shape)
        #print(x.T.shape)
        dXdP = np.dot(X.T, dNdP)
        dXdQ = np.dot(X.T, dNdQ)
        dXdR = np.dot(X.T, dNdR)
        dYdP = np.dot(Y.T, dNdP)
        dYdQ = np.dot(Y.T, dNdQ)
        dYdR = np.dot(Y.T, dNdR)
        dZdP = np.dot(Z.T, dNdP)
        dZdQ = np.dot(Z.T, dNdQ)
        dZdR = np.dot(Z.T, dNdR)
        #print(dXdP.shape)
        dxdP = np.dot(x.T, dNdP)
        dxdQ = np.dot(x.T, dNdQ)
        dxdR = np.dot(x.T, dNdR)
        dydP = np.dot(y.T, dNdP)
        dydQ = np.dot(y.T, dNdQ)
        dydR = np.dot(y.T, dNdR)
        dzdP = np.dot(z.T, dNdP)
        dzdQ = np.dot(z.T, dNdQ)
        dzdR = np.dot(z.T, dNdR)

        dudP = np.dot(u.T, dNdP)
        dudQ = np.dot(u.T, dNdQ)
        dudR = np.dot(u.T, dNdR)
        dvdP = np.dot(v.T, dNdP)
        dvdQ = np.dot(v.T, dNdQ)
        dvdR = np.dot(v.T, dNdR)
        dwdP = np.dot(w.T, dNdP)
        dwdQ = np.dot(w.T, dNdQ)
        dwdR = np.dot(w.T, dNdR)

        FXP = np.array([np.concatenate([dXdP, dXdQ, dXdR]), np.concatenate([dYdP, dYdQ, dYdR]), np.concatenate([dZdP, dZdQ, dZdR])])[:,:,0]
        FxP = np.array([np.concatenate([dxdP, dxdQ, dxdR]), np.concatenate([dydP, dydQ, dydR]), np.concatenate([dzdP, dzdQ, dzdR])])[:,:,0]
        FuP = np.array([np.concatenate([dudP, dudQ, dudR]), np.concatenate([dvdP, dvdQ, dvdR]), np.concatenate([dwdP, dwdQ, dwdR])])[:,:,0]
        #print(FXP)
        #print(FXP.shape)
        FxX = np.linalg.inv(FXP).dot(FxP)
        FuX = FuP.dot(np.linalg.inv(FXP))
        I = np.eye(20)
        XFXP = np.concatenate([np.concatenate([dXdP * I, dYdP * I, dZdP * I],axis=1),
                         np.concatenate([dXdQ * I, dYdQ * I, dZdQ * I],axis=1),
                        np.concatenate([dXdR * I, dYdR * I, dZdR * I],axis=1)])
        #print(XFXP.shape)
        vec = np.linalg.inv(XFXP).dot(np.array([dNdP, dNdQ, dNdR]).reshape(-1,1))
        dNdX = vec[0:20].reshape(-1,1)
        dNdY = vec[20:40].reshape(-1,1)
        dNdZ = vec[40:60].reshape(-1,1)

        dudX = FuX[0, 0]
        dudY = FuX[0, 1]
        dudZ = FuX[0, 2]
        dvdX = FuX[1, 0]
        dvdY = FuX[1, 1]
        dvdZ = FuX[1, 2]
        dwdX = FuX[2, 0]
        dwdY = FuX[2, 1]
        dwdZ = FuX[2, 2]
        #print(dNdX.shape)
        ZERO = np.zeros([20,1])

        BL = np.concatenate([np.concatenate([dNdX, ZERO, ZERO, dNdY, ZERO, dNdZ],axis=1),
                       np.concatenate([ZERO, dNdY, ZERO, dNdX, dNdZ, ZERO],axis=1),
                        np.concatenate([ZERO, ZERO, dNdZ, ZERO, dNdY, dNdX],axis=1)])
        BNN = np.concatenate([np.concatenate([dudX*dNdX, dudY*dNdY, dudZ*dNdZ],axis=1),
                       np.concatenate([dvdX*dNdX, dvdY*dNdY, dvdZ*dNdZ],axis=1),
                        np.concatenate([dwdX*dNdX, dwdY*dNdY, dwdZ*dNdZ],axis=1)])
        BNS = np.concatenate([np.concatenate([dudX*dNdY+dudY*dNdX, dudY*dNdZ+dudZ*dNdY, dudZ*dNdX+dudX*dNdZ],axis=1),
                       np.concatenate([dvdX*dNdY+dvdY*dNdX, dvdY*dNdZ+dvdZ*dNdY, dvdZ*dNdX+dvdX*dNdZ],axis=1),
                        np.concatenate([dwdX*dNdY+dwdY*dNdX, dwdY*dNdZ+dwdZ*dNdY, dwdZ*dNdX+dwdX*dNdZ],axis=1)])
        #print(BL.shape)
        #print(np.concatenate([BNN, BNS],axis=1).shape)
        B = BL + np.concatenate([BNN, BNS],axis=1)
        J = np.linalg.det(FXP)
        F = np.copy(FxX)
        return B,F,J

    @staticmethod
    def mk_dNdP_vec(P,Q,R,DIM):
        V = np.zeros([20,1])
        if(DIM is 1):
            V[0] = 1 / 8 * (1 - Q) * (1 - R) * (1 + 2 * P + Q + R)
            V[1] = -1 / 4 * (1 - R ** 2) * (1 - Q)
            V[2] = 1 / 8 * (1 - Q) * (1 + R) * (1 + 2 * P + Q - R)
            V[3] = -1 / 4 * (1 - Q ** 2) * (1 - R)
            V[4] = -1 / 4 * (1 - Q ** 2) * (1 + R)
            V[4] = 1 / 8 * (1 + Q) * (1 - R) * (1 + 2 * P - Q + R)
            V[5] = -1 / 4 * (1 - R ** 2) * (1 + Q)
            V[7] = 1 / 8 * (1 + Q) * (1 + R) * (1 + 2 * P - Q - R)
            V[8] = -1 / 2 * P * (1 - Q) * (1 - R)
            V[9] = -1 / 2 * P * (1 - Q) * (1 + R)
            V[10] = -1 / 2 * P * (1 + Q) * (1 - R)
            V[11] = -1 / 2 * P * (1 + Q) * (1 + R)
            V[12] = -1 / 8 * (1 - Q) * (1 - R) * (1 - 2 * P + Q + R)
            V[13] = 1 / 4 * (1 - R ** 2) * (1 - Q)
            V[14] = -1 / 8 * (1 - Q) * (1 + R) * (1 - 2 * P + Q - R)
            V[15] = 1 / 4 * (1 - Q ** 2) * (1 - R)
            V[16] = 1 / 4 * (1 - Q ** 2) * (1 + R)
            V[17] = -1 / 8 * (1 + Q) * (1 - R) * (1 - 2 * P - Q + R)
            V[18] = 1 / 4 * (1 - R ** 2) * (1 + Q)
            V[19] = -1 / 8 * (1 + Q) * (1 + R) * (1 - 2 * P - Q - R)
        elif(DIM is 2):
            V[0] = 1 / 8 * (1 - P) * (1 - R) * (1 + P + 2 * Q + R)
            V[1] = -1 / 4 * (1 - R ** 2) * (1 - P)
            V[2] = 1 / 8 * (1 - P) * (1 + R) * (1 + P + 2 * Q - R)
            V[3] = -1 / 2 * Q * (1 - P) * (1 - R)
            V[4] = -1 / 2 * Q * (1 - P) * (1 + R)
            V[5] = -1 / 8 * (1 - P) * (1 - R) * (1 + P - 2 * Q + R)
            V[6] = 1 / 4 * (1 - R ** 2) * (1 - P)
            V[7] = -1 / 8 * (1 - P) * (1 + R) * (1 + P - 2 * Q - R)
            V[8] = -1 / 4 * (1 - P ** 2) * (1 - R)
            V[9] = -1 / 4 * (1 - P ** 2) * (1 + R)
            V[10] = 1 / 4 * (1 - P ** 2) * (1 - R)
            V[11] = 1 / 4 * (1 - P ** 2) * (1 + R)
            V[12] = 1 / 8 * (1 + P) * (1 - R) * (1 - P + 2 * Q + R)
            V[13] = -1 / 4 * (1 - R ** 2) * (1 + P)
            V[14] = 1 / 8 * (1 + P) * (1 + R) * (1 - P + 2 * Q - R)
            V[15] = -1 / 2 * Q * (1 + P) * (1 - R)
            V[16] = -1 / 2 * Q * (1 + P) * (1 + R)
            V[17] = -1 / 8 * (1 + P) * (1 - R) * (1 - P - 2 * Q + R)
            V[18] = 1 / 4 * (1 - R ** 2) * (1 + P)
            V[19] = -1 / 8 * (1 + P) * (1 + R) * (1 - P - 2 * Q - R)
        elif(DIM is 3):
            V[0] = 1 / 8 * (1 - P) * (1 - Q) * (1 + P + Q + 2 * R)
            V[1] = -1 / 2 * R * (1 - P) * (1 - Q)
            V[2] = -1 / 8 * (1 - P) * (1 - Q) * (1 + P + Q - 2 * R)
            V[3] = -1 / 4 * (1 - Q ** 2) * (1 - P)
            V[4] = 1 / 4 * (1 - Q ** 2) * (1 - P)
            V[5] = 1 / 8 * (1 - P) * (1 + Q) * (1 + P - Q + 2 * R)
            V[6] = -1 / 2 * R * (1 - P) * (1 + Q)
            V[7] = -1 / 8 * (1 - P) * (1 + Q) * (1 + P - Q - 2 * R)
            V[8] = -1 / 4 * (1 - P ** 2) * (1 - Q)
            V[9] = 1 / 4 * (1 - P ** 2) * (1 - Q)
            V[10] = -1 / 4 * (1 - P ** 2) * (1 + Q)
            V[11] = 1 / 4 * (1 - P ** 2) * (1 + Q)
            V[12] = 1 / 8 * (1 + P) * (1 - Q) * (1 - P + Q + 2 * R)
            V[13] = -1 / 2 * R * (1 + P) * (1 - Q)
            V[14] = -1 / 8 * (1 + P) * (1 - Q) * (1 - P + Q - 2 * R)
            V[15] = -1 / 4 * (1 - Q ** 2) * (1 + P)
            V[16] = 1 / 4 * (1 - Q ** 2) * (1 + P)
            V[17] = 1 / 8 * (1 + P) * (1 + Q) * (1 - P - Q + 2 * R)
            V[18] = -1 / 2 * R * (1 + P) * (1 + Q)
            V[19] = -1 / 8 * (1 + P) * (1 + Q) * (1 - P - Q - 2 * R)
        return V

    @staticmethod
    def mk_Mises_stress(Sij,F):
        Sij = Sij.reshape(-1,1)
        #print(Sij.shape)
        S = np.array([np.concatenate([Sij[0], Sij[3], Sij[5]]),
                       np.concatenate([Sij[3], Sij[1], Sij[4]]),
                        np.concatenate([Sij[5], Sij[4], Sij[2]])])
        #print(S.shape)
        K = np.dot(np.dot(F, S), F.T)
        SP = K - np.sum(np.diag(K)) * np.eye(3) / 3
        J2 = 1 / 2 * np.sum(np.diag(SP) ** 2) + SP[0, 1] ** 2 + SP[1, 2] ** 2 + SP[2, 0] ** 2
        ms = np.sqrt(3 * J2)
        return ms

    @staticmethod
    def mk_global_mtr(K,Udof,dUb,Fdof,dFb):
        #print(dFb.shape)
        #print(dUb.shape)
        #print(Fdof.shape)
        #print(Udof.shape)
        rhs = dFb[:, 0].reshape(-1,1) - np.dot(K[np.ix_(np.ravel(Fdof)-1, Udof-1)], dUb[:, 0].reshape(-1,1))
        Kb = K[np.ix_(np.ravel(Fdof)-1, np.ravel(Fdof)-1)]
        return Kb, rhs

    @staticmethod
    def lgwt(N,a,b):
        N = N - 1
        N1 = N + 1
        N2 = N + 2
        xu = np.linspace(-1,1,N1).T
        #print("test")
        #print(np.cos((2*np.array(range(0,N+1)).T+1)*np.pi).shape)
        #print(((0.27/N1)*np.sin(np.pi*xu*N/N2)).shape)
        y = np.cos((2*np.array(range(0,N+1)).T+1)*np.pi/(2*N+2)) + (0.27/N1)*np.sin(np.pi*xu*N/N2)
        y = y.reshape(-1,1)
        #print(y.shape)
        #print(y)
        L = np.zeros([N1,N2])
        Lp = np.zeros([N1,N2])
        y0 = 2
        while np.max(np.abs(y-y0)) > sys.float_info.epsilon:
            L[:, 0] = 1
            L[:, 1] = y[:,0]
            #print(L)
            for k in range(2,N1+1):
                L[:, k] = (((2 * k - 1) * y * L[:, k-1].reshape(-1,1) - (k - 1) * L[:, k - 2].reshape(-1,1)) / k)[:,0]

            Lp = N2 * (L[:, N1-1].reshape(-1,1) - y * L[:, N2-1].reshape(-1,1)) / (1 - y ** 2)
            Lp = Lp.reshape(Lp.shape[0],1)
            y0 = np.copy(y.reshape(-1,1))
            y = y0 - L[:, N2-1].reshape(-1,1) / Lp
            #print(L.shape)
            #print(y.shape)
        p = (a * (1 - y) + b * (1 + y)) / 2
        w = (b - a) / ((1 - y ** 2) * Lp ** 2) * (N2 / N1) ** 2
        return p,w

