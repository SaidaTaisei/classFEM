import param

if __name__=="__main__":
    param = param.param()
    node = []
    for x in range(2*param.x_divide_num):
        for y in range(2 * param.x_divide_num):
            for z in range(2 * param.x_divide_num):
                if(x%2 + y%2 + z%2 < 2):
                    pos_x = x * param.x_length / 2 / param.x_divide_num
                    pos_y = y * param.y_length / 2 / param.y_divide_num
                    pos_z = z * param.z_length / 2 / param.z_divide_num
                    node.append([x, y, z, pos_x, pos_y, pos_z, 0, 0, 0])
    elem = []
    for x in range(2*param.x_divide_num):
        for y in range(2 * param.x_divide_num):
            for z in range(2 * param.x_divide_num):
                x1 = (2 * x - 1) - 1
                x2 = (2 * x - 1)
                x3 = (2 * x - 1) + 1
                y1 = (2 * y - 1) - 1
                y2 = (2 * y - 1)
                y3 = (2 * y - 1) + 1
                z1 = (2 * z - 1) - 1
                z2 = (2 * z - 1)
                z3 = (2 * z - 1) + 1