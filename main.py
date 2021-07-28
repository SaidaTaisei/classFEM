import model

if __name__=="__main__":
    model = model.model()
    model.create_node()
    model.create_elem()
    model.set_bcs()