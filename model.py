import numpy as np

class model:
    def __init__(self):
        self.x_length = 10.0
        self.y_length = 1.0
        self.z_length = 1.0
        self.x_elemnum = 20
        self.y_elemnum = 2
        self.z_elemnum = 2
        self.E = 200000
        self.nu = 0.2
        self.H = 3000
        self.fy = 5000
        self.step_num = 20

    def create_node(self):
        self.node = []
        for x in range(2*self.x_elemnum+1):
            for y in range(2*self.y_elemnum+1):
                for z in range(2*self.z_elemnum+1):
                    if(x%2 + y%2 + z%2 < 2):
                        x_pos = x * self.x_length / 2.0 / self.x_elemnum
                        y_pos = y * self.y_length / 2.0 / self.y_elemnum
                        z_pos = z * self.z_length / 2.0 / self.z_elemnum
                        self.node.append([x,y,z,x_pos,y_pos,z_pos,0,0,0])
                        print(x_pos)
        self.node = np.array(self.node)
        print(self.node)
        print(self.node.shape)

    def create_elem(self):
        self.elem = []
        for x in range(1,self.x_elemnum+1):
            for y in range(1,self.y_elemnum+1):
                for z in range(1,self.z_elemnum+1):
                    x1 = (2 * x - 1) - 1
                    x2 = (2 * x - 1)
                    x3 = (2 * x - 1) + 1
                    y1 = (2 * y - 1) - 1
                    y2 = (2 * y - 1)
                    y3 = (2 * y - 1) + 1
                    z1 = (2 * z - 1) - 1
                    z2 = (2 * z - 1)
                    z3 = (2 * z - 1) + 1

                    index = []
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y1, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y1, z2], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y1, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y2, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y2, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y3, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y3, z2], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x1, y3, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x2, y1, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x2, y1, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x2, y3, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x2, y3, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y1, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y1, z2], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y1, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y2, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y2, z3], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y3, z1], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y3, z2], axis=1) == 3)[0] + 1)
                    index.append(np.where(np.count_nonzero(self.node[:, 0:3] == [x3, y3, z3], axis=1) == 3)[0] + 1)
                    index = np.array(index)
                    self.elem.append(index.T)
        self.elem = np.array(self.elem)
        self.elem = self.elem.reshape(self.elem.shape[0],self.elem.shape[-1])
        print(self.elem)
        print(self.elem.shape)

    def set_bcs(self):
        self.nnum = self.node.shape[0]
        self.enum = self.elem.shape[0]
        self.dof = 3 * self.nnum
        self.bc_node_0 = np.where((self.node[:, 0] == 0) & (self.node[:, 2] == 0))[0] + 1
        bc_dof_0 = np.concatenate([self.bc_node_0,self.nnum+self.bc_node_0,2*self.nnum+self.bc_node_0])
        dub_0 = np.zeros([bc_dof_0.shape[0],self.step_num])
        self.bc_node_1 = np.where((self.node[:, 0] == 2*self.x_elemnum) & (self.node[:, 2] == 0))[0] + 1
        bc_dof_1 = np.concatenate([self.bc_node_1, 2*self.nnum + self.bc_node_1])
        dub_1 = np.zeros([bc_dof_1.shape[0], self.step_num])
        self.bc_node_2 = np.where((self.node[:, 0] == self.x_elemnum) & (self.node[:, 2] == 2*self.z_elemnum))[0] + 1
        bc_dof_2 = 2*self.nnum+self.bc_node_2
        self.inc = -0.2*np.ones([1,self.step_num])
        dub_2 = -0.2*np.ones([bc_dof_2.shape[0],self.inc.shape[1]])
        self.udof = np.concatenate([bc_dof_0,bc_dof_1,bc_dof_2])
        self.dub = np.concatenate([dub_0,dub_1,dub_2])
        self.fdof = np.array(range(self.dof))+1
        self.fdof = self.fdof.reshape([self.fdof.shape[0],1])
        self.fdof = np.delete(self.fdof,self.udof)
        self.fdof = self.fdof.reshape([self.fdof.shape[0],1])
        self.dfb = np.zeros(self.fdof.shape[0],self.step_num)
        print(self.fdof.shape)
        print(self.bc_dof_1)
        print(self.dub_2.shape)