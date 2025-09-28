
class PDES:
    eqs = []
    ivars = []
    sp_vars = []
    funcs = []
    ic = []
    n_sp = 1
    n_temp = 1
    def __init__(self, pdes, sp_vars, n_sp=1, n_temp=1):
        self.eqs = [pde.eq for pde in pdes]
        self.ivars = list([pde.ivar[i] for pde in pdes for i in range(len(pde.ivar))])
        self.funcs = list([pde.func[i] for pde in pdes for i in range(len(pde.func))])
        self.ic = list([pde.ic[i] for pde in pdes for i in range(len(pde.ic))])
        self.sp_vars = sp_vars
        self.n_sp = n_sp
        self.n_temp = n_temp
        pass
    
    
    def xs(self,vars):
        nvars = vars.copy()
        for i in range(len(nvars)): nvars[i] = f'XX{i}'
        return nvars

