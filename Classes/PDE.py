import sympy as sp
import numpy as np
from FuncAux import symbol_references
class PDE:
    
    eq = ''
    func = ''
    expr_ic = ''
    ivar = []
    disc_n = []
    ivar_boundary = []
    ic = []
    
    def __init__(self, eq, func, ivar, disc_n, ivar_boundary,expr_ic):
        self.eq = eq
        self.func = func     
        self.expr_ic = expr_ic
        self.ivar = ivar     
        self.disc_n = disc_n
        self.ivar_boundary = ivar_boundary
        self.ic = self.ic_calc()
    
    def ic_calc(self):
        ivar_symbols = symbol_references(self.ivar)
        expr = sp.parse_expr(self.expr_ic)

        # Gera grade de pontos para cada variável
        grids = []
        for i, var in enumerate(self.ivar):
            a, b = self.ivar_boundary[i]
            n = self.disc_n[i]
            grids.append(np.linspace(a, b, n))

        
        # Cria malha multidimensional (como meshgrid)
        mesh = np.meshgrid(*grids, indexing='ij')
        
        # Converte expressão simbólica para função NumPy
        f_ic = sp.lambdify(self.ivar, expr, modules='numpy')
        
        # Avalia a função na malha de pontos
        ic = f_ic(*mesh)
        
        return ic.flatten().tolist()  # array NumPy com os ic nos pontos da malha


