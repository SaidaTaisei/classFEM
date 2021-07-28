import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

class visualizer:
    @staticmethod
    def visualize_matplotlib(model):
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot(model.node[:,3], model.node[:,4], model.node[:,5], "k.")
        print(model.node.shape[0])
        print(model.elem.shape[0])
        for ee in range(1,model.elem.shape[0]+1):
            side = []
            side.append(np.array([model.node[model.elem[ee-1,0]-1,3:6], model.node[model.elem[ee-1,1]-1,3:6], model.node[model.elem[ee-1,2]-1,3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 0]-1, 3:6], model.node[model.elem[ee - 1, 3]-1, 3:6],
                                  model.node[model.elem[ee - 1, 5]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 2]-1, 3:6], model.node[model.elem[ee - 1, 4]-1, 3:6],
                                  model.node[model.elem[ee - 1, 7]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 5]-1, 3:6], model.node[model.elem[ee - 1, 6]-1, 3:6],
                                  model.node[model.elem[ee - 1, 7]-1, 3:6]]))

            side.append(np.array([model.node[model.elem[ee - 1, 0]-1, 3:6], model.node[model.elem[ee - 1, 8]-1, 3:6],
                                  model.node[model.elem[ee - 1, 12]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 2]-1, 3:6], model.node[model.elem[ee - 1, 9]-1, 3:6],
                                  model.node[model.elem[ee - 1, 14]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 5]-1, 3:6], model.node[model.elem[ee - 1, 10]-1, 3:6],
                                  model.node[model.elem[ee - 1, 17]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 7]-1, 3:6], model.node[model.elem[ee - 1, 11]-1, 3:6],
                                  model.node[model.elem[ee - 1, 19]-1, 3:6]]))

            side.append(np.array([model.node[model.elem[ee - 1, 12]-1, 3:6], model.node[model.elem[ee - 1, 13]-1, 3:6],
                                  model.node[model.elem[ee - 1, 14]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 12]-1, 3:6], model.node[model.elem[ee - 1, 15]-1, 3:6],
                                  model.node[model.elem[ee - 1, 17]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 14]-1, 3:6], model.node[model.elem[ee - 1, 16]-1, 3:6],
                                  model.node[model.elem[ee - 1, 19]-1, 3:6]]))
            side.append(np.array([model.node[model.elem[ee - 1, 17]-1, 3:6], model.node[model.elem[ee - 1, 18]-1, 3:6],
                                  model.node[model.elem[ee - 1, 19]-1, 3:6]]))
            for ii in range(1,13):
                ax.plot(side[ii-1][:,0], side[ii-1][:,1], side[ii-1][:,2], "b-")
        for ee in range(1,model.elem.shape[0]+1):
            side = []
            side.append(np.array([model.node[model.elem[ee-1,0]-1,6:9], model.node[model.elem[ee-1,1]-1,6:9], model.node[model.elem[ee-1,2]-1,6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 0]-1, 6:9], model.node[model.elem[ee - 1, 3]-1, 6:9],
                                  model.node[model.elem[ee - 1, 5]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 2]-1, 6:9], model.node[model.elem[ee - 1, 4]-1, 6:9],
                                  model.node[model.elem[ee - 1, 7]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 5]-1, 6:9], model.node[model.elem[ee - 1, 6]-1, 6:9],
                                  model.node[model.elem[ee - 1, 7]-1, 6:9]]))

            side.append(np.array([model.node[model.elem[ee - 1, 0]-1, 6:9], model.node[model.elem[ee - 1, 8]-1, 6:9],
                                  model.node[model.elem[ee - 1, 12]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 2]-1, 6:9], model.node[model.elem[ee - 1, 9]-1, 6:9],
                                  model.node[model.elem[ee - 1, 14]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 5]-1, 6:9], model.node[model.elem[ee - 1, 10]-1, 6:9],
                                  model.node[model.elem[ee - 1, 17]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 7]-1, 6:9], model.node[model.elem[ee - 1, 11]-1, 6:9],
                                  model.node[model.elem[ee - 1, 19]-1, 6:9]]))

            side.append(np.array([model.node[model.elem[ee - 1, 12]-1, 6:9], model.node[model.elem[ee - 1, 13]-1, 6:9],
                                  model.node[model.elem[ee - 1, 14]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 12]-1, 6:9], model.node[model.elem[ee - 1, 15]-1, 6:9],
                                  model.node[model.elem[ee - 1, 17]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 14]-1, 6:9], model.node[model.elem[ee - 1, 16]-1, 6:9],
                                  model.node[model.elem[ee - 1, 19]-1, 6:9]]))
            side.append(np.array([model.node[model.elem[ee - 1, 17]-1, 6:9], model.node[model.elem[ee - 1, 18]-1, 6:9],
                                  model.node[model.elem[ee - 1, 19]-1, 6:9]]))
            for ii in range(1,13):
                ax.plot(side[ii-1][:,0], side[ii-1][:,1], side[ii-1][:,2], "r-")
        #rmax = np.max()
        plt.show()