import model
import numpy as np
import mk
import solver
import visualizer
import pandas as pd

if __name__=="__main__":
    model = model.model()
    model.create_node()
    model.create_elem()
    model.set_bcs()

    max_itr = 10
    lgn = 2
    max_err = 10**-6

    bc_num, step_num = model.dub.shape
    nnum = model.nnum
    enum = model.enum
    dof = model.dof
    U = np.zeros([dof,step_num])
    dU = np.zeros([dof,step_num])
    dU[model.udof-1,:] = model.dub

    x = np.concatenate([model.node[:, 3].reshape(model.node.shape[0], 1), np.zeros([nnum, step_num - 1])], axis=1)
    y = np.concatenate([model.node[:, 4].reshape(model.node.shape[0], 1), np.zeros([nnum, step_num - 1])], axis=1)
    z = np.concatenate([model.node[:, 5].reshape(model.node.shape[0], 1), np.zeros([nnum, step_num - 1])], axis=1)
    E = model.E * np.ones([enum, 1])
    nu = model.nu * np.ones([enum, 1])
    H = model.H * np.ones([enum, 1])
    fy = model.fy * np.ones([enum, 1])
    S = np.zeros([enum*lgn**3, 6])

    smax = []
    smax = np.array(smax)
    for ee in range(1,enum+1):
        if (ee is 1):
            smax = fy[ee-1,0]*np.ones([lgn**3,1])
        else:
            smax = np.concatenate([smax,fy[ee-1,0]*np.ones([lgn**3,1])])
    print(smax.shape)
    for tt in range(1,step_num):
        xt = x[:, tt - 1].reshape(-1,1)
        yt = y[:, tt - 1].reshape(-1,1)
        zt = z[:, tt - 1].reshape(-1,1)
        Ut = U[:, tt - 1].reshape(-1,1)
        dUt = dU[:, tt - 1].reshape(-1,1)
        dUbt = model.dub[:, tt - 1].reshape(-1,1)
        dFbt = model.dfb[:, tt - 1].reshape(-1,1)
        dS = np.zeros([enum*lgn**3,6])
        sMax = 0.0
        for itr in range(1,max_itr+1):
            print(str(tt) + " - " + str(itr))
            dUt0 = np.copy(dUt)
            K,dS,sMax = mk.mk.mk_stif_mtr(model.node,model.elem,xt,yt,zt,Ut,dUt,E,H,nu,S+dS,smax)
            dUt[model.fdof-1,0], FLG = solver.solver.solve_dU(K, model.udof, dUbt, model.fdof, dFbt)
            ERR = np.linalg.norm(dUt-dUt0)/np.linalg.norm(dUt0)
            print("ERR : "+str(ERR))
            if(ERR <= max_err):
                break
        dU[:, tt-1] = dUt[:,0]
        U[:, tt] = U[:, tt-1]+dU[:, tt-1]
        x[:, tt] = model.node[:, 3] + U[0: nnum, tt]
        y[:, tt] = model.node[:, 4] + U[nnum+np.array(range(0, nnum)), tt]
        z[:, tt] = model.node[:, 5] + U[2*nnum+np.array(range(0, nnum)), tt]
        smax = np.copy(sMax)
        S = S + dS

    u = U[0:nnum,:]
    v = U[nnum:2 * nnum,:]
    w = U[2 * nnum:3 * nnum,:]
    model.node[:, 6] = model.node[:, 3] + u[:, -2]
    model.node[:, 7] = model.node[:, 4] + v[:, -2]
    model.node[:, 8] = model.node[:, 5] + w[:, -2]
    node_data = pd.DataFrame(model.node)
    node_data.to_csv("result_node.csv")
    visualizer.visualizer.visualize_matplotlib(model)
